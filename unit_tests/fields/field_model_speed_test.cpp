/*
 * field_model_speed_test.cpp
 *
 *  Created on: Sept 21, 2020
 *      Author: David Flanderka
 *
 *  Tests evaluation of FieldModel
 */

#define TEST_USE_MPI
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>
#include <mesh_constructor.hh>

#include "fields/eval_points.hh"
#include "fields/eval_subset.hh"
#include "fields/field_value_cache.hh"
#include "fields/field_values.hh"
#include "fields/field_set.hh"
#include "tools/unit_si.hh"
#include "fields/field.hh"
#include "fields/field_model.hh"
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

class Balance;


/****************************************************************************************
 *                 Speed test of FieldModel evaluation
 *
 * Results:
 * Mesh 27936 elements, 100 assemblation loops
 * Checked GenericAssembly with active bulk integral only vs. with all active integrals
 * (times in brackets are )
 *
 *                           bulk                     all
 * add_integrals_to_patch   18.99                    42.38
 * create_patch              2.97                    18.30
 * cache_update             10.30  11.63* -8.43**    31.38  34.19* -26.38**
 * *) times of cache_update with functions passed to model ( not struct and operator() )
 * **) times of FieldConstant fields of field_evaluate_const_test
 *
 ****************************************************************************************/

#ifdef FLOW123D_RUN_UNIT_BENCHMARKS

static const unsigned int profiler_loop = 100;

using Sclr = double;
using Vect = arma::vec3;
using Tens = arma::mat33;

// Functors of FieldModels
struct fn_model_scalar {
    inline Sclr operator() (Vect vec) {
        return arma::norm(vec, 2);
    }
};

struct fn_model_vector {
	inline Vect operator() (Sclr scal, Vect vec) {
        return scal * vec + vec;
    }
};

struct fn_model_tensor {
	inline Tens operator() (Sclr scal, Tens ten) {
        return 0.5 * scal * ten;
    }
};


class FieldModelSpeedTest : public testing::Test {
public:
    class EqData : public FieldSet {
    public:
        EqData() : order(2) {
            *this += vector_const_field
                        .name("vector_const_field")
                        .description("Constant vector.")
                        .input_default("0.0")
                        .flags_add(in_main_matrix)
                        .units( UnitSI().kg(3).m() );
            *this += scalar_const_field
                        .name("scalar_const_field")
                        .description("Constant scalar")
                        .units( UnitSI().m() );
            *this += tensor_const_field
                        .name("tensor_const_field")
                        .description("Constant tensor")
                        .units( UnitSI::dimensionless() )
                        .flags_add(in_main_matrix);
            *this += vector_field
                        .name("vector_field")
                        .description("Velocity vector.")
                        .input_default("0.0")
                        .flags_add(in_main_matrix)
                        .units( UnitSI().kg(3).m() );
            *this += scalar_field
                        .name("scalar_field")
                        .description("Pressure head")
                        .input_default("0.0")
                        .units( UnitSI().m() );
            *this += tensor_field
                        .name("tensor_field")
                        .description("")
                        .input_default("0.0")
                        .units( UnitSI::dimensionless() )
                        .flags_add(in_main_matrix);

        }

        /// Initialiye FieldModels
        void initialize() {
            scalar_field.set(Model<3, FieldValue<3>::Scalar>::create(fn_model_scalar(), vector_const_field), 0.0);
            vector_field.set(Model<3, FieldValue<3>::VectorFixed>::create(fn_model_vector(), scalar_const_field, vector_const_field), 0.0);
            tensor_field.set(Model<3, FieldValue<3>::TensorFixed>::create(fn_model_tensor(), scalar_const_field, tensor_const_field), 0.0);
        }

    	/// Polynomial order of finite elements.
    	unsigned int order;

    	// constant fields, we need these fields to create models
        Field<3, FieldValue<3>::Scalar > scalar_const_field;
        Field<3, FieldValue<3>::VectorFixed > vector_const_field;
        Field<3, FieldValue<3>::TensorFixed > tensor_const_field;

        // computing fields of FieldModels
        Field<3, FieldValue<3>::Scalar > scalar_field;
        Field<3, FieldValue<3>::VectorFixed > vector_field;
        Field<3, FieldValue<3>::TensorFixed > tensor_field;
    };

    FieldModelSpeedTest() : tg_(0.0, 1.0) {
    	Profiler::instance();

        FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
        Profiler::instance();
        PetscInitialize(0,PETSC_NULL,PETSC_NULL,PETSC_NULL);

        data_ = std::make_shared<EqData>();
        mesh_ = mesh_full_constructor("{ mesh_file=\"mesh/test_27936_elem.msh\", optimize_mesh=false }");
        dh_ = std::make_shared<DOFHandlerMultiDim>(*mesh_);
    }

    ~FieldModelSpeedTest() {
        Profiler::uninitialize();
    }

    static Input::Type::Record & get_input_type() {
        return IT::Record("SomeEquation","")
                .declare_key("data", IT::Array(
                        IT::Record("SomeEquation_Data", FieldCommon::field_descriptor_record_description("SomeEquation_Data") )
                        .copy_keys( FieldModelSpeedTest::EqData().make_field_descriptor_type("SomeEquation") )
                        .declare_key("scalar_const_field", FieldAlgorithmBase< 3, FieldValue<3>::Scalar >::get_input_type_instance(), "" )
                        .declare_key("vector_const_field", FieldAlgorithmBase< 3, FieldValue<3>::VectorFixed >::get_input_type_instance(), "" )
                        .declare_key("tensor_const_field", FieldAlgorithmBase< 3, FieldValue<3>::TensorFixed >::get_input_type_instance(), "" )
//                        .declare_key("scalar_field", FieldAlgorithmBase< 3, FieldValue<3>::Scalar >::get_input_type_instance(), "" )
//                        .declare_key("vector_field", FieldAlgorithmBase< 3, FieldValue<3>::VectorFixed >::get_input_type_instance(), "" )
//                        .declare_key("tensor_field", FieldAlgorithmBase< 3, FieldValue<3>::TensorFixed >::get_input_type_instance(), "" )
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
        data_->initialize();
        data_->set_time(tg_.step(), LimitSide::right);
    }


	void profiler_output() {
		static ofstream os( FilePath("speed_eval_model_test.log", FilePath::output_file) );
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
    typedef typename FieldModelSpeedTest::EqData EqFields;
    typedef typename FieldModelSpeedTest::EqData EqData;

    static constexpr const char * name() { return "AssemblyModelTest"; }

    /// Constructor.
    AssemblyDimTest(EqFields *eq_fields, EqData *eq_data, PatchFEValues<3> *fe_values)
    : AssemblyBase<dim>(eq_data->order, fe_values), eq_fields_(eq_fields), eq_data_(eq_data) {
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
    scalar_const_field: !FieldConstant
      value: 0.5
    vector_const_field: [1, 2, 3]
    tensor_const_field: [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
)YAML";


TEST_F(FieldModelSpeedTest, speed_test) {
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
