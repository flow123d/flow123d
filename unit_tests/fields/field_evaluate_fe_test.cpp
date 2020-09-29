/*
 * field_evaluate_fe_test.cpp
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
#include "mesh/mesh.h"
#include "mesh/accessors.hh"
#include "input/input_type.hh"
#include "input/accessors.hh"
#include "input/reader_to_storage.hh"
#include "system/sys_profiler.hh"
#include "coupling/generic_assembly.hh"
#include "coupling/assembly_base.hh"


class FieldEvalFETest : public testing::Test {

public:
    class EqData : public FieldSet, public ElementCacheMap {
    public:
        EqData() {
            *this += vector_field
                        .name("vector_field")
                        .description("Velocity vector.")
                        .input_default("0.0")
                        .flags_add(in_main_matrix)
                        .units( UnitSI().kg(3).m() );
            *this += scalar_field
                        .name("scalar_field")
                        .description("Pressure head")
                        .units( UnitSI().m() );
            *this += tensor_field
                        .name("tensor_field")
                        .description("")
                        .units( UnitSI::dimensionless() )
                        .flags_add(in_main_matrix);

            // Asumme following types:
            eval_points_ = std::make_shared<EvalPoints>();
            Quadrature *q_bulk = new QGauss(3, 2);
            Quadrature *q_side = new QGauss(2, 2);
            mass_eval = eval_points_->add_bulk<3>(*q_bulk );
            side_eval = eval_points_->add_edge<3>(*q_side );
            // ngh_side_eval = ...
            this->init(eval_points_);
        }

        void register_eval_points() {
            unsigned int reg_idx = computed_dh_cell_.elm().region_idx().idx();
            for (auto p : mass_eval->points(computed_dh_cell_, this) ) {
                EvalPointData epd(reg_idx, computed_dh_cell_.elm_idx(), p.eval_point_idx());
                this->eval_point_data_.push_back(epd);
            }

            for (DHCellSide cell_side : computed_dh_cell_.side_range()) {
            	for( DHCellSide edge_side : cell_side.edge_sides() ) {
                    unsigned int reg_idx = edge_side.element().region_idx().idx();
                    for (auto p : side_eval->points(edge_side, this) ) {
                        EvalPointData epd(reg_idx, edge_side.elem_idx(), p.eval_point_idx());
                        this->eval_point_data_.push_back(epd);
                    }
                }
            }
            this->eval_point_data_.make_permanent();
        }

        void update_cache() {
            this->register_eval_points();
            this->create_patch();
            this->cache_update(*this);
            this->finish_elements_update();
        }

        void realocate_cache() {
            this->cache_reallocate(*this);
        }


        // fields
        Field<3, FieldValue<3>::Scalar > scalar_field;
        Field<3, FieldValue<3>::VectorFixed > vector_field;
        Field<3, FieldValue<3>::TensorFixed > tensor_field;
        std::shared_ptr<EvalPoints> eval_points_;
        std::shared_ptr<BulkIntegral> mass_eval;
        std::shared_ptr<EdgeIntegral> side_eval;
        //std::shared_ptr<CouplingIntegral> ngh_side_eval;
        DHCellAccessor computed_dh_cell_;
    };

    FieldEvalFETest() {
        FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
        Profiler::instance();
        PetscInitialize(0,PETSC_NULL,PETSC_NULL,PETSC_NULL);

        data_ = std::make_shared<EqData>();
        mesh_ = mesh_full_constructor("{mesh_file=\"mesh/cube_2x1.msh\"}");
        dh_ = std::make_shared<DOFHandlerMultiDim>(*mesh_);
    }

    ~FieldEvalFETest() {}

    static Input::Type::Record & get_input_type() {
        return IT::Record("SomeEquation","")
                .declare_key("data", IT::Array(
                        IT::Record("SomeEquation_Data", FieldCommon::field_descriptor_record_description("SomeEquation_Data") )
                        .copy_keys( FieldEvalFETest::EqData().make_field_descriptor_type("SomeEquation") )
                        .declare_key("scalar_field", FieldAlgorithmBase< 3, FieldValue<3>::Scalar >::get_input_type_instance(), "" )
                        .declare_key("vector_field", FieldAlgorithmBase< 3, FieldValue<3>::VectorFixed >::get_input_type_instance(), "" )
                        .declare_key("tensor_field", FieldAlgorithmBase< 3, FieldValue<3>::TensorFixed >::get_input_type_instance(), "" )
                        .close()
                        ), IT::Default::obligatory(), ""  )
                .close();
    }

    void read_input(const string &input) {
        // read input string
        Input::ReaderToStorage reader( input, get_input_type(), Input::FileFormat::format_YAML );
        Input::Record in_rec=reader.get_root_interface<Input::Record>();

        TimeGovernor tg(0.0, 1.0);

        //data.set_components(component_names);        // set number of substances posibly read from elsewhere

        static std::vector<Input::Array> inputs;
        unsigned int input_last = inputs.size(); // position of new item
        inputs.push_back( in_rec.val<Input::Array>("data") );

        data_->set_mesh(*mesh_);
        data_->set_input_list( inputs[input_last], tg );
        data_->set_time(tg.step(), LimitSide::right);
        data_->realocate_cache();
    }


    std::shared_ptr<EqData> data_;
    Mesh * mesh_;
    std::shared_ptr<DOFHandlerMultiDim> dh_;
};


TEST_F(FieldEvalFETest, evaluate) {
    string eq_data_input = R"YAML(
    data:
      - region: 3D left
        time: 0.0
        scalar_field: !FieldFE
          mesh_data_file: mesh/cube_2x1.msh
          field_name: scalar
        vector_field: !FieldFE
          mesh_data_file: mesh/cube_2x1.msh
          field_name: vector
        tensor_field: !FieldFE
          mesh_data_file: mesh/cube_2x1.msh
          field_name: tensor
      - region: 3D right
        time: 0.0
        scalar_field: !FieldFE
          mesh_data_file: mesh/cube_2x1.msh
          field_name: scalar
        vector_field: !FieldFE
          mesh_data_file: mesh/cube_2x1.msh
          field_name: vector
        tensor_field: !FieldFE
          mesh_data_file: mesh/cube_2x1.msh
          field_name: tensor
    )YAML";
	this->read_input(eq_data_input);

    std::vector<unsigned int> cell_idx = {3, 4, 5, 9};
    std::vector<double>       expected_scalar = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.1, 1.2, 1.3};
    std::vector<arma::vec3>   expected_vector = {{1, 4, 7}, {2, 5, 8}, {3, 6, 9}, {1, 4, 7}};
    std::vector<arma::mat33>  expected_tensor = {{3.1, 3.4, 3.7, 3.2, 3.5, 3.8, 3.3, 3.6, 3.9}, {4.1, 4.4, 4.7, 4.2, 4.5, 4.8, 4.3, 4.6, 4.9},
                                                 {5.1, 5.4, 5.7, 5.2, 5.5, 5.8, 5.3, 5.6, 5.9}, {9.1, 9.4, 9.7, 9.2, 9.5, 9.8, 9.3, 9.6, 9.9}};
    for (unsigned int i=0; i<cell_idx.size(); ++i) {
    	data_->start_elements_update();
    	data_->computed_dh_cell_ = DHCellAccessor(dh_.get(), cell_idx[i]);  // element ids stored to cache: (3 -> 2,3,4), (4 -> 3,4,5,10), (5 -> 0,4,5,11), (10 -> 8,9,10)
        data_->update_cache();

        // Bulk integral, no sides, no permutations.
        for(BulkPoint q_point: data_->mass_eval->points(data_->computed_dh_cell_, data_.get())) {
            EXPECT_EQ(expected_scalar[data_->computed_dh_cell_.elm_idx()], data_->scalar_field(q_point));
            EXPECT_ARMA_EQ(expected_vector[i], data_->vector_field(q_point));
            EXPECT_ARMA_EQ(expected_tensor[i], data_->tensor_field(q_point));
            /* // Extracting the cached values.
            double cs = cross_section(q_point);

            // Following would be nice to have. Not clear how to
            // deal with more then single element as fe_values have its own cache that has to be updated.
            auto base_fn_grad = presssure_field_fe.base_value(q_point);
            loc_matrix += outer_product((cs * base_fn_grad),  base_fn_grad) */
        }

        // Side integrals.
        // FieldFE<..> conc;
        for (DHCellSide side : data_->computed_dh_cell_.side_range()) {
        	for(DHCellSide el_ngh_side : side.edge_sides()) {
           	    // vector of local side quadrature points in the correct side permutation
        	    Range<EdgePoint> side_points = data_->side_eval->points(side, data_.get());
        	    for (EdgePoint side_p : side_points) {
                    EXPECT_EQ(expected_scalar[data_->computed_dh_cell_.elm_idx()], data_->scalar_field(side_p));
                    EXPECT_ARMA_EQ(expected_vector[i], data_->vector_field(side_p));
                    EXPECT_ARMA_EQ(expected_tensor[i], data_->tensor_field(side_p));
                    EdgePoint ngh_p = side_p.point_on(el_ngh_side);
                    EXPECT_EQ(expected_scalar[el_ngh_side.cell().elm_idx()], data_->scalar_field(ngh_p));
        	        //loc_mat += cross_section(side_p) * sigma(side_p) *
        		    //    (conc.base_value(side_p) * velocity(side_p)
        		    //    + conc.base_value(ngh_p) * velocity(ngh_p)) * side_p.normal() / 2;
                }
            }
        }
    }

}


/****************************************************************************************
 *                 Speed test of FieldFE evaluation
 *
 * Results:
 * Mesh 27936 elements, 50 assemblation loops
 * Checked GenericAssembly with active bulk integral only vs. with all active integrals
 *
 *                           bulk            all
 * add_integrals_to_patch   19.10 (19.28)   44.12  (44.38)
 * create_patch              3.14  (3.12)   19.28  (19.32)
 * cache_update             58.30 (41.16)  314.10 (267.20)
 * (times are multiplied by 2 for simple compare with other speed tests with 100
 * assemblation loops - FieldConstant and FieldModel)
 * (values in brackets are for inlined functions get_loc_dof_indices and shape_value_component)
 *
 ****************************************************************************************/

#ifdef FLOW123D_RUN_UNIT_BENCHMARKS

static const unsigned int profiler_loop = 50;


class FieldFESpeedTest : public testing::Test {
public:
    class EqData : public FieldSet {
    public:
        EqData() : order(2) {
            *this += vector_field
                        .name("vector_field")
                        .description("Velocity vector.")
                        .input_default("0.0")
                        .flags_add(in_main_matrix)
                        .units( UnitSI().kg(3).m() );
            *this += scalar_field
                        .name("scalar_field")
                        .description("Pressure head")
                        .units( UnitSI().m() );
            *this += tensor_field
                        .name("tensor_field")
                        .description("")
                        .units( UnitSI::dimensionless() )
                        .flags_add(in_main_matrix);

        }

    	/// Polynomial order of finite elements.
    	unsigned int order;

    	// fields
        Field<3, FieldValue<3>::Scalar > scalar_field;
        Field<3, FieldValue<3>::VectorFixed > vector_field;
        Field<3, FieldValue<3>::TensorFixed > tensor_field;
    };

    FieldFESpeedTest() : tg_(0.0, 1.0) {
    	Profiler::instance();

        FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
        Profiler::instance();
        PetscInitialize(0,PETSC_NULL,PETSC_NULL,PETSC_NULL);

        data_ = std::make_shared<EqData>();
        mesh_ = mesh_full_constructor("{mesh_file=\"mesh/test_27936_elem.msh\"}");
        dh_ = std::make_shared<DOFHandlerMultiDim>(*mesh_);
    }

    ~FieldFESpeedTest() {
        Profiler::uninitialize();
    }

    static Input::Type::Record & get_input_type() {
        return IT::Record("SomeEquation","")
                .declare_key("data", IT::Array(
                        IT::Record("SomeEquation_Data", FieldCommon::field_descriptor_record_description("SomeEquation_Data") )
                        .copy_keys( FieldFESpeedTest::EqData().make_field_descriptor_type("SomeEquation") )
                        .declare_key("scalar_field", FieldAlgorithmBase< 3, FieldValue<3>::Scalar >::get_input_type_instance(), "" )
                        .declare_key("vector_field", FieldAlgorithmBase< 3, FieldValue<3>::VectorFixed >::get_input_type_instance(), "" )
                        .declare_key("tensor_field", FieldAlgorithmBase< 3, FieldValue<3>::TensorFixed >::get_input_type_instance(), "" )
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

    	//auto scalar_vec = data_->scalar_field.get_field_fe()->vec();
    	//auto vector_vec = data_->vector_field.get_field_fe()->vec();
    	//auto tensor_vec = data_->tensor_field.get_field_fe()->vec();
    	//for (unsigned int i=0; i<scalar_vec.size(); ++i) {
    	//    scalar_vec[i] = 1 + i % 9;
    	//}
    	//for (unsigned int i=0; i<vector_vec.size(); ++i) {
    	//	vector_vec[i] = 1 + i % 9;
    	//}
    	//for (unsigned int i=0; i<tensor_vec.size(); ++i) {
    	//	tensor_vec[i] = 1 + i % 9;
    	//}

        data_->set_time(tg_.step(), LimitSide::right);
    }


	void profiler_output() {
		static ofstream os( FilePath("speed_eval_fe_test.log", FilePath::output_file) );
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
    typedef typename FieldFESpeedTest::EqData EqDataDG;

    /// Constructor.
    AssemblyDimTest(EqDataDG *data)
    : AssemblyBase<dim>(data->order), data_(data) {}

    void initialize(FMT_UNUSED std::shared_ptr<Balance> balance) {
    }

    void reallocate_cache(const ElementCacheMap &cache_map) override
    {
        data_->cache_reallocate(cache_map);
    }

    /// Data object shared with Test class
    EqDataDG *data_;
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


TEST_F(FieldFESpeedTest, speed_test) {
	this->read_input(eq_data_input_speed);

	std::shared_ptr<Balance> balance;
	GenericAssembly< AssemblyDimTest > ga_bulk(data_.get(), balance, ActiveIntegrals::bulk);
	START_TIMER("assemble_bulk");
	for (unsigned int i=0; i<profiler_loop; ++i)
		ga_bulk.assemble(this->dh_, this->tg_.step());
	END_TIMER("assemble_bulk");

	GenericAssembly< AssemblyDimTest > ga_all(data_.get(), balance,
	        (ActiveIntegrals::bulk | ActiveIntegrals::edge | ActiveIntegrals::coupling | ActiveIntegrals::boundary) );
	START_TIMER("assemble_all_integrals");
	for (unsigned int i=0; i<profiler_loop; ++i)
		ga_all.assemble(this->dh_, this->tg_.step());
	END_TIMER("assemble_all_integrals");

	this->profiler_output();
}


#endif // FLOW123D_RUN_UNIT_BENCHMARKS
