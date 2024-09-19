/*
 * field_const_speed_test.cpp
 *
 *  Created on: Dec 03, 2019
 *      Author: David Flanderka
 *
 *  Speed tests of FieldConstant
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
#include "tools/unit_si.hh"
#include "quadrature/quadrature.hh"
#include "quadrature/quadrature_lib.hh"
#include "fem/dofhandler.hh"
#include "fem/dh_cell_accessor.hh"
#include "mesh/mesh.h"
#include "mesh/accessors.hh"
#include "input/input_type.hh"
#include "input/accessors.hh"
#include "input/reader_to_storage.hh"
#include "system/sys_profiler.hh"
#include "coupling/generic_assembly.hh"
#include "coupling/assembly_base.hh"
#include "coupling/equation.hh"

class Balance;


/****************************************************************************************
 *                 Speed test of FieldConstant evaluation
 *
 * Results:
 * Mesh 27936 elements, 100 assemblation loops
 * Checked GenericAssembly with active bulk integral only vs. with all active integrals
 *
 *                           bulk      all
 * add_integrals_to_patch   19.12    43.95
 * create_patch              3.00    18.59
 * cache_update              8.43    26.38
 *
 ****************************************************************************************/

#ifdef FLOW123D_RUN_UNIT_BENCHMARKS

static const unsigned int profiler_loop = 100;

// allow running tests of all type of fields (unset allow runnig only vector field)
#define ALL_FIELDS


class FieldConstantSpeedTest : public testing::Test {
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

    FieldConstantSpeedTest() : tg_(0.0, 1.0) {
    	Profiler::instance();

        FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
        Profiler::instance();
        PetscInitialize(0,PETSC_NULL,PETSC_NULL,PETSC_NULL);

        data_ = std::make_shared<EqData>();
        data_->add_coords_field();
        mesh_ = mesh_full_constructor("{ mesh_file=\"mesh/test_27936_elem.msh\", optimize_mesh=false }");
        dh_ = std::make_shared<DOFHandlerMultiDim>(*mesh_);
    }

    ~FieldConstantSpeedTest() {
        Profiler::uninitialize();
    }

    static Input::Type::Record & get_input_type() {
        return IT::Record("SomeEquation","")
                .declare_key("data", IT::Array(
                        IT::Record("SomeEquation_Data", FieldCommon::field_descriptor_record_description("SomeEquation_Data") )
                        .copy_keys( FieldEvalConstantTest::EqData().make_field_descriptor_type("SomeEquation") )
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


	void profiler_output() {
		static ofstream os( FilePath("speed_eval_const_test.log", FilePath::output_file) );
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
    typedef typename FieldConstantSpeedTest::EqData EqFields;
    typedef typename FieldConstantSpeedTest::EqData EqData;

    static constexpr const char * name() { return "AssemblyConstantTest"; }

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
    scalar_field: !FieldConstant
      value: 0.5
    vector_field: [1, 2, 3]
    tensor_field: [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
)YAML";


TEST_F(FieldConstantSpeedTest, speed_test) {
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

	this->profiler_output();
}


#endif // FLOW123D_RUN_UNIT_BENCHMARKS
