/*
 * field_time_function_test.cpp
 *
 *  Created on: Dec 18, 2012
 *      Author: jb
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
#include "fields/field_constant.hh"
#include "fields/field_time_function.hh"
#include "fields/table_function.hh"
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

class Balance;


class FieldTimeFunctionTest : public testing::Test {

public:
    class EqData : public FieldSet, public ElementCacheMap {
    public:
        EqData() {
            *this += vector_field
                        .name("vector_field")
                        .description("vector field.")
                        .input_default("0.0")
                        .flags_add(in_main_matrix)
                        .units( UnitSI().kg(3).m() );
            *this += scalar_field
                        .name("scalar_field")
                        .description("scalar field")
                        .units( UnitSI().m() );
            *this += tensor_field
                        .name("tensor_field")
                        .description("Tensor field")
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
            for (auto p : mass_eval->points(this->position_in_cache(computed_dh_cell_.elm_idx()), this) ) {
                this->eval_point_data_.emplace_back(reg_idx, computed_dh_cell_.elm_idx(), p.eval_point_idx(), computed_dh_cell_.local_idx());
            }

            for (DHCellSide cell_side : computed_dh_cell_.side_range()) {
            	for( DHCellSide edge_side : cell_side.edge_sides() ) {
                    unsigned int reg_idx = edge_side.element().region_idx().idx();
                    for (auto p : side_eval->points(edge_side, this) ) {
                        this->eval_point_data_.emplace_back(reg_idx, edge_side.elem_idx(), p.eval_point_idx(), edge_side.cell().local_idx());
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

    FieldTimeFunctionTest()
    : tg(0.0, 0.5) {
        FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
        Profiler::instance();
        PetscInitialize(0,PETSC_NULL,PETSC_NULL,PETSC_NULL);

        data_ = std::make_shared<EqData>();
        data_->add_coords_field();
        mesh_ = mesh_full_constructor("{ mesh_file=\"mesh/cube_2x1.msh\", optimize_mesh=false }");
        dh_ = std::make_shared<DOFHandlerMultiDim>(*mesh_);
    }

    ~FieldTimeFunctionTest() {}

    static Input::Type::Record & get_input_type() {
        return IT::Record("SomeEquation","")
                .declare_key("data", IT::Array(
                        IT::Record("SomeEquation_Data", FieldCommon::field_descriptor_record_description("SomeEquation_Data") )
                        .copy_keys( FieldTimeFunctionTest::EqData().make_field_descriptor_type("SomeEquation") )
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
        data_->set_input_list( inputs[input_last], tg );
        data_->set_time(tg.step(), LimitSide::right);
        data_->cache_reallocate( *(data_.get()), *(data_.get()) );
    }


    std::shared_ptr<EqData> data_;
    Mesh * mesh_;
    std::shared_ptr<DOFHandlerMultiDim> dh_;
    TimeGovernor tg;
};

string eq_data_input = R"YAML(
data:
  - region: BULK
    time: 0.0
    scalar_field: !FieldTimeFunction
      time_function:
        - - 0.0
          - 0.5
        - - 1.0
          - 1.0
        - - 2.0
          - 3.0
    vector_field: !FieldTimeFunction
      time_function:
        - - 0.0
          - [0.5, 1.5, 2.0]
        - - 1.0
          - [1.0, 1.5, 1.0]
        - - 2.0
          - [3.0, 1.5, 5.0]
    tensor_field: !FieldTimeFunction
      time_function:
        - - 0.0
          - [ [1,3,4], [0,3,4], [1,6,6] ]
        - - 1.0
          - [ [3,3,6], [0,2,5], [2,7,4] ]
        - - 2.0
          - [ [5,3,4], [2,3,7], [5,6,3] ]
)YAML";

TEST_F(FieldTimeFunctionTest, evaluate) {
    this->read_input(eq_data_input);

    unsigned int             n_times = 6;   // time steps: { 0.0, 0.5, 1.0, 1.5, 2.0, 2.5 }
    unsigned int             elem_idx = 4;  // test only on one element
    std::vector<double>      expected_scalar = {0.5, 0.75, 1.0, 2.0, 3.0, 3.0};
    std::vector<arma::vec3>  expected_vector = {{0.5, 1.5, 2.0}, {0.75, 1.5, 1.5}, {1.0, 1.5, 1.0}, {2.0, 1.5, 3.0}, {3.0, 1.5, 5.0}, {3.0, 1.5, 5.0}};
    std::vector<arma::mat33> expected_tensor = {{1.0, 0.0, 1.0, 3.0, 3.0, 6.0, 4.0, 4.0, 6.0}, {2.0, 0.0, 1.5, 3.0, 2.5, 6.5, 5.0, 4.5, 5.0},
                                                {3.0, 0.0, 2.0, 3.0, 2.0, 7.0, 6.0, 5.0, 4.0}, {4.0, 1.0, 3.5, 3.0, 2.5, 6.5, 5.0, 6.0, 3.5},
                                                {5.0, 2.0, 5.0, 3.0, 3.0, 6.0, 4.0, 7.0, 3.0}, {5.0, 2.0, 5.0, 3.0, 3.0, 6.0, 4.0, 7.0, 3.0}};

    for (unsigned int i=0; i<n_times; ++i) {
    	data_->set_time(tg.step(), LimitSide::right);
    	data_->start_elements_update();
    	data_->computed_dh_cell_ = DHCellAccessor(dh_.get(), elem_idx);  // element ids stored to cache: (4 -> 3,4,5,10)
        data_->update_cache();

        // Bulk integral, no sides.
        for( BulkPoint q_point: data_->mass_eval->points(data_->position_in_cache(data_->computed_dh_cell_.elm_idx()), data_.get()) ) {
            EXPECT_EQ(expected_scalar[i], data_->scalar_field(q_point));
            EXPECT_ARMA_EQ(expected_vector[i], data_->vector_field(q_point));
            EXPECT_ARMA_EQ(expected_tensor[i], data_->tensor_field(q_point));
        }

        // Side integrals.
        for (DHCellSide side : data_->computed_dh_cell_.side_range()) {
        	for(DHCellSide el_ngh_side : side.edge_sides()) {
           	    // vector of local side quadrature points
        	    Range<EdgePoint> side_points = data_->side_eval->points(side, data_.get());
        	    for (EdgePoint side_p : side_points) {
                    EXPECT_EQ(expected_scalar[i], data_->scalar_field(side_p));
                    EXPECT_ARMA_EQ(expected_vector[i], data_->vector_field(side_p));
                    EXPECT_ARMA_EQ(expected_tensor[i], data_->tensor_field(side_p));
                    EdgePoint ngh_p = side_p.point_on(el_ngh_side);
                    EXPECT_EQ(expected_scalar[i], data_->scalar_field(ngh_p));
                }
            }
        }

        tg.next_time();
    }

}
