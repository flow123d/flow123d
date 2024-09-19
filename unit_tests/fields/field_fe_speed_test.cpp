/*
 * field_fe_speed_test.cpp
 *
 *  Created on: Apr 07, 2020
 *      Author: David Flanderka
 *
 *  Tests evaluation of FieldFE
 */

#define TEST_USE_MPI
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>
#include <mesh_constructor.hh>
#include "arma_expect.hh"

#include "fields/eval_points.hh"
#include "fields/eval_subset.hh"
#include "fields/field_value_cache.hh"
#include "fields/field_values.hh"
#include "fields/field_set.hh"
#include "fields/field_fe.hh"
#include "tools/unit_si.hh"
#include "fields/bc_field.hh"
#include "quadrature/quadrature.hh"
#include "quadrature/quadrature_lib.hh"
#include "fem/dofhandler.hh"
#include "fem/dh_cell_accessor.hh"
#include "fem/fe_rt.hh"
#include "mesh/mesh.h"
#include "mesh/accessors.hh"
#include "input/input_type.hh"
#include "input/accessors.hh"
#include "input/reader_to_storage.hh"
#include "system/sys_profiler.hh"
#include "coupling/generic_assembly.hh"
#include "coupling/assembly_base.hh"
#include "coupling/equation.hh"


/****************************************************************************************
 *                 Speed test of FieldFE evaluation
 *
 * Results:
 * Mesh 27936 elements, 50 assemblation loops
 * Checked GenericAssembly with active bulk integral only vs. with all active integrals
 *
 *                        |     scalar + vector + tensor      |        only vector field          |          FieldConstant            |
 *                        |   FieldFE file  |  FieldFE FE RT  |   FieldFE file  |  FieldFE FE RT  | scalar+vec+tens | only vec field  |
 *                        |   bulk |   all  |   bulk |   all  |   bulk |   all  |   bulk |   all  |   bulk |   all  |   bulk |   all  |
 * add_integrals_to_patch |  19.16 |  44.88 |  19.44 |  44.68 |  18.96 |  43.10 |  19.62 |  45.08 |  19.03 |  45.35 |  19.23 |  43.54 |
 * create_patch           |   3.20 |  19.52 |   3.12 |  19.48 |   3.02 |  18.48 |   3.06 |  19.48 |   3.01 |  20.06 |   3.05 |  18.38 |
 * cache_update           |  86.96 | 398.80 |  87.16 | 403.68 |  24.08 | 106.28 |  24.66 | 109.22 |   7.91 |  25.01 |   7.57 |  23.36 |
 * (times are multiplied by 2 for simple compare with other speed tests with 100
 * assemblation loops - FieldConstant and FieldModel)
 * (values in brackets are for inlined functions get_loc_dof_indices and shape_value_component)
 *
 ****************************************************************************************/

#ifdef FLOW123D_RUN_UNIT_BENCHMARKS

static const unsigned int profiler_loop = 50;

// allow running tests of all type of fields (unset allow runnig only vector field)
#define ALL_FIELDS

class FieldFESpeedTest : public testing::Test {
public:
    class EqData : public FieldSet, public EqDataBase {
    public:
        EqData() : EqDataBase(2) {
            *this += vector_field
                        .name("vector_field")
                        .description("Velocity vector.")
                        .input_default("0.0")
                        .flags_add(in_main_matrix)
                        .units( UnitSI().kg(3).m() );
#ifdef ALL_FIELDS
            *this += scalar_field
                        .name("scalar_field")
                        .description("Pressure head")
                        .units( UnitSI().m() );
            *this += tensor_field
                        .name("tensor_field")
                        .description("")
                        .units( UnitSI::dimensionless() )
                        .flags_add(in_main_matrix);
#endif // ALL_FIELDS
        }

    	// fields
        Field<3, FieldValue<3>::VectorFixed > vector_field;
#ifdef ALL_FIELDS
    	Field<3, FieldValue<3>::Scalar > scalar_field;
        Field<3, FieldValue<3>::TensorFixed > tensor_field;
#endif // ALL_FIELDS
    };

    FieldFESpeedTest() : tg_(0.0, 1.0) {
    	Profiler::instance();

        FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
        Profiler::instance();
        PetscInitialize(0,PETSC_NULL,PETSC_NULL,PETSC_NULL);

        data_ = std::make_shared<EqData>();
        mesh_ = mesh_full_constructor("{ mesh_file=\"mesh/test_27936_elem.msh\", optimize_mesh=false }");
    	MixedPtr<FE_RT0> fe_rt0;
        std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>(mesh_, fe_rt0);
        dh_ = std::make_shared<DOFHandlerMultiDim>(*mesh_);
        dh_->distribute_dofs(ds);
    }

    ~FieldFESpeedTest() {
        Profiler::uninitialize();
    }

    static Input::Type::Record & get_input_type() {
        return IT::Record("SomeEquation","")
                .declare_key("data", IT::Array(
                        IT::Record("SomeEquation_Data", FieldCommon::field_descriptor_record_description("SomeEquation_Data") )
                        .copy_keys( FieldFESpeedTest::EqData().make_field_descriptor_type("SomeEquation") )
                        .declare_key("vector_field", FieldAlgorithmBase< 3, FieldValue<3>::VectorFixed >::get_input_type_instance(), "" )
#ifdef ALL_FIELDS
                        .declare_key("scalar_field", FieldAlgorithmBase< 3, FieldValue<3>::Scalar >::get_input_type_instance(), "" )
                        .declare_key("tensor_field", FieldAlgorithmBase< 3, FieldValue<3>::TensorFixed >::get_input_type_instance(), "" )
#endif // ALL_FIELDS
                        .close()
                        ), IT::Default::obligatory(), ""  )
                .close();
    }

    void read_input(const string &input) {
        // read input string
        Input::ReaderToStorage reader( input, get_input_type(), Input::FileFormat::format_YAML );
        Input::Record in_rec=reader.get_root_interface<Input::Record>();

        //data.set_components(component_names);        // set number of substances posibly read from elsewhere

        static std::vector<Input::Array> inputs;
        unsigned int input_last = inputs.size(); // position of new item
        inputs.push_back( in_rec.val<Input::Array>("data") );

        data_->set_mesh(*mesh_);
        data_->set_input_list( inputs[input_last], tg_ );

        data_->set_time(tg_.step(), LimitSide::right);
    }

    void create_fe_fields() {
    	VectorMPI * vector_vec = new VectorMPI(dh_->distr()->lsize() * 3);
    	auto vector_ptr = create_field_fe<3, FieldValue<3>::VectorFixed>(dh_, vector_vec);
        data_->vector_field.set(vector_ptr, 0.0);
    	for (unsigned int i=0; i<vector_vec->size(); ++i) {
    		vector_vec->set( i, (i % 10 + 0.5) );
    	}

#ifdef ALL_FIELDS
        VectorMPI * scalar_vec = new VectorMPI(dh_->distr()->lsize());
        auto scalar_ptr = create_field_fe<3, FieldValue<3>::Scalar>(dh_, scalar_vec);
        data_->scalar_field.set(scalar_ptr, 0.0);
    	for (unsigned int i=0; i<scalar_vec->size(); ++i) {
    	    scalar_vec->set( i, (1 + i % 9) );
    	}

    	VectorMPI * tensor_vec = new VectorMPI(dh_->distr()->lsize() * 9);
        auto tensor_ptr = create_field_fe<3, FieldValue<3>::TensorFixed>(dh_, tensor_vec);
        data_->tensor_field.set(tensor_ptr, 0.0);
    	for (unsigned int i=0; i<tensor_vec->size(); ++i) {
    		tensor_vec->set( i, (1 + i % 98) * 0.1 );
    	}
#endif // ALL_FIELDS

        data_->set_mesh(*mesh_);
        data_->set_time(tg_.step(), LimitSide::right);
    }


	void profiler_output(std::string field_fill) {
		static ofstream os( FilePath("speed_eval_fe_" + field_fill + "_test.log", FilePath::output_file) );
		Profiler::instance()->output(MPI_COMM_WORLD, os);
		os << "" << std::setfill('=') << setw(80) << "" << std::setfill(' ') << endl << endl;
	}

    std::shared_ptr<EqData> data_;
    Mesh * mesh_;
    std::shared_ptr<DOFHandlerMultiDim> dh_;
    TimeGovernor tg_;
};

template <unsigned int dim>
class AssemblyDimTest : public AssemblyBase<dim> {
public:
    typedef typename FieldFESpeedTest::EqData EqFields;
    typedef typename FieldFESpeedTest::EqData EqData;

    static constexpr const char * name() { return "AssemblyFETest"; }

    /// Constructor.
    AssemblyDimTest(EqFields *eq_fields, EqData *eq_data)
    : AssemblyBase<dim>(eq_data->quad_order), eq_fields_(eq_fields), eq_data_(eq_data) {
        this->active_integrals_ = (ActiveIntegrals::bulk | ActiveIntegrals::edge | ActiveIntegrals::coupling);
        this->used_fields_.set_mesh( *eq_fields_->mesh() );
        this->used_fields_ += *eq_fields_;
    }

    void initialize(ElementCacheMap *element_cache_map) {
        this->element_cache_map_ = element_cache_map;
    }

    /// Data object shared with Test class
    EqFields *eq_fields_;
    EqData *eq_data_;

    /// Sub field set contains fields used in calculation.
    FieldSet used_fields_;
};

string eq_data_input_speed = R"YAML(
data:
  - region: ALL
    time: 0.0
    scalar_field: !FieldFE
      mesh_data_file: mesh/test_27936_elem.msh
      field_name: scalar
      default_value: 0.0
    vector_field: !FieldFE
      mesh_data_file: mesh/test_27936_elem.msh
      field_name: vector
      default_value: 0.1
    tensor_field: !FieldFE
      mesh_data_file: mesh/test_27936_elem.msh
      field_name: tensor
      default_value: 0.2
)YAML";


TEST_F(FieldFESpeedTest, read_from_input_test) {
	this->read_input(eq_data_input_speed);

	GenericAssembly< AssemblyDimTest > ga_bulk(data_.get(), data_.get());
	START_TIMER("assemble_bulk");
	for (unsigned int i=0; i<profiler_loop; ++i)
		ga_bulk.assemble(this->dh_);
	END_TIMER("assemble_bulk");

	GenericAssembly< AssemblyDimTest > ga_all(data_.get(), data_.get());
	START_TIMER("assemble_all_integrals");
	for (unsigned int i=0; i<profiler_loop; ++i)
		ga_all.assemble(this->dh_);
	END_TIMER("assemble_all_integrals");

	this->profiler_output("file");
}

TEST_F(FieldFESpeedTest, rt_field_fe_test) {
	this->read_input(eq_data_input_speed);

	GenericAssembly< AssemblyDimTest > ga_bulk(data_.get(), data_.get());
	START_TIMER("assemble_bulk");
	for (unsigned int i=0; i<profiler_loop; ++i)
		ga_bulk.assemble(this->dh_);
	END_TIMER("assemble_bulk");

	GenericAssembly< AssemblyDimTest > ga_all(data_.get(), data_.get());
	START_TIMER("assemble_all_integrals");
	for (unsigned int i=0; i<profiler_loop; ++i)
		ga_all.assemble(this->dh_);
	END_TIMER("assemble_all_integrals");

	this->profiler_output("rt");
}


#endif // FLOW123D_RUN_UNIT_BENCHMARKS
