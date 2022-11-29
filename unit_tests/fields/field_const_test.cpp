/*
 * field_const_test.cpp
 *
 *  Created on: Dec 03, 2019
 *      Author: David Flanderka
 *
 *  Tests evaluation of FieldConstant
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

class Balance;


class FieldEvalConstantTest : public testing::Test {

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
            *this += tensor_1
                        .name("tensor_1")
                        .description("Tensor defined by single value.")
                        .units( UnitSI::dimensionless() )
                        .flags_add(in_main_matrix);
            *this += tensor_2
                        .name("tensor_2")
                        .description("Tensor defined by diagonal values.")
                        .units( UnitSI::dimensionless() )
                        .flags_add(in_main_matrix);
            *this += tensor_3
                        .name("tensor_3")
                        .description("Tensor defined by 6 values as symmetric matrix.")
                        .units( UnitSI::dimensionless() )
                        .flags_add(in_main_matrix);
            *this += tensor_4
                        .name("tensor_4")
                        .description("Tensor defined by 9 values as full matrix.")
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
        Field<3, FieldValue<3>::TensorFixed > tensor_1;
        Field<3, FieldValue<3>::TensorFixed > tensor_2;
        Field<3, FieldValue<3>::TensorFixed > tensor_3;
        Field<3, FieldValue<3>::TensorFixed > tensor_4;
        std::shared_ptr<EvalPoints> eval_points_;
        std::shared_ptr<BulkIntegral> mass_eval;
        std::shared_ptr<EdgeIntegral> side_eval;
        //std::shared_ptr<CouplingIntegral> ngh_side_eval;
        DHCellAccessor computed_dh_cell_;
    };

    FieldEvalConstantTest() {
        FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
        Profiler::instance();
        PetscInitialize(0,PETSC_NULL,PETSC_NULL,PETSC_NULL);

        data_ = std::make_shared<EqData>();
        data_->add_coords_field();
        mesh_ = mesh_full_constructor("{ mesh_file=\"mesh/cube_2x1.msh\", optimize_mesh=false }");
        dh_ = std::make_shared<DOFHandlerMultiDim>(*mesh_);
    }

    ~FieldEvalConstantTest() {}

    static Input::Type::Record & get_input_type() {
        return IT::Record("SomeEquation","")
                .declare_key("data", IT::Array(
                        IT::Record("SomeEquation_Data", FieldCommon::field_descriptor_record_description("SomeEquation_Data") )
                        .copy_keys( FieldEvalConstantTest::EqData().make_field_descriptor_type("SomeEquation") )
                        .declare_key("scalar_field", FieldAlgorithmBase< 3, FieldValue<3>::Scalar >::get_input_type_instance(), "" )
                        .declare_key("vector_field", FieldAlgorithmBase< 3, FieldValue<3>::VectorFixed >::get_input_type_instance(), "" )
                        .declare_key("tensor_1", FieldAlgorithmBase< 3, FieldValue<3>::TensorFixed >::get_input_type_instance(), "" )
                        .declare_key("tensor_2", FieldAlgorithmBase< 3, FieldValue<3>::TensorFixed >::get_input_type_instance(), "" )
                        .declare_key("tensor_3", FieldAlgorithmBase< 3, FieldValue<3>::TensorFixed >::get_input_type_instance(), "" )
                        .declare_key("tensor_4", FieldAlgorithmBase< 3, FieldValue<3>::TensorFixed >::get_input_type_instance(), "" )
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
        data_->cache_reallocate( *(data_.get()), *(data_.get()) );
    }


    std::shared_ptr<EqData> data_;
    Mesh * mesh_;
    std::shared_ptr<DOFHandlerMultiDim> dh_;
};

string eq_data_input = R"YAML(
data:
  - region: 3D left
    time: 0.0
    scalar_field: !FieldConstant
      value: 0.5
    vector_field: [1, 2, 3]
    tensor_1: 3.14
    tensor_2: [1, 2, 3]
    tensor_3: [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
    tensor_4: [ [0,1,2], [3,4,5], [6,7,8] ]
  - region: 3D right
    time: 0.0
    scalar_field: !FieldConstant
      value: 15
      unit: "dm"
    vector_field: [4, 5, 6]
    tensor_1: 2.72
    tensor_2: [1.2, 2.3, 3.4]
    tensor_3: [2.1, 2.2, 2.3, 2.4, 2.5, 2.6]
    tensor_4: [ [1,2,3], [4,5,6], [7,8,9] ]
)YAML";

TEST_F(FieldEvalConstantTest, evaluate) {
    this->read_input(eq_data_input);

    std::vector<unsigned int> cell_idx = {3, 4, 5, 9};
    std::vector<double>       expected_scalar = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5};
    std::vector<arma::vec3>   expected_vector = {{1, 2, 3}, {1, 2, 3}, {1, 2, 3}, {4, 5, 6}};
    std::vector<arma::mat33>  expected_tens_1 = {{3.14, 0, 0, 0, 3.14, 0, 0, 0, 3.14}, {3.14, 0, 0, 0, 3.14, 0, 0, 0, 3.14},
                                                 {3.14, 0, 0, 0, 3.14, 0, 0, 0, 3.14}, {2.72, 0, 0, 0, 2.72, 0, 0, 0, 2.72}};
    std::vector<arma::mat33>  expected_tens_2 = {{1.0, 0, 0, 0, 2.0, 0, 0, 0, 3.0}, {1.0, 0, 0, 0, 2.0, 0, 0, 0, 3.0},
                                                 {1.0, 0, 0, 0, 2.0, 0, 0, 0, 3.0}, {1.2, 0, 0, 0, 2.3, 0, 0, 0, 3.4}};
    std::vector<arma::mat33>  expected_tens_3 = {{0.1, 0.2, 0.3, 0.2, 0.4, 0.5, 0.3, 0.5, 0.6}, {0.1, 0.2, 0.3, 0.2, 0.4, 0.5, 0.3, 0.5, 0.6},
                                                 {0.1, 0.2, 0.3, 0.2, 0.4, 0.5, 0.3, 0.5, 0.6}, {2.1, 2.2, 2.3, 2.2, 2.4, 2.5, 2.3, 2.5, 2.6}};
    std::vector<arma::mat33>  expected_tens_4 = {{0.0, 3, 6, 1, 4, 7, 2, 5, 8}, {0.0, 3, 6, 1, 4, 7, 2, 5, 8},
                                                 {0.0, 3, 6, 1, 4, 7, 2, 5, 8}, {1.0, 4, 7, 2, 5, 8, 3, 6, 9}};

    for (unsigned int i=0; i<cell_idx.size(); ++i) {
    	data_->start_elements_update();
    	data_->computed_dh_cell_ = DHCellAccessor(dh_.get(), cell_idx[i]);  // element ids stored to cache: (3 -> 2,3,4), (4 -> 3,4,5,10), (5 -> 0,4,5,11), (10 -> 8,9,10)
        data_->update_cache();

        // Bulk integral, no sides.
        for( BulkPoint q_point: data_->mass_eval->points(data_->position_in_cache(data_->computed_dh_cell_.elm_idx()), data_.get()) ) {
            EXPECT_EQ(expected_scalar[data_->computed_dh_cell_.elm_idx()], data_->scalar_field(q_point));
            EXPECT_ARMA_EQ(expected_vector[i], data_->vector_field(q_point));
            EXPECT_ARMA_EQ(expected_tens_1[i], data_->tensor_1(q_point));
            EXPECT_ARMA_EQ(expected_tens_2[i], data_->tensor_2(q_point));
            EXPECT_ARMA_EQ(expected_tens_3[i], data_->tensor_3(q_point));
            EXPECT_ARMA_EQ(expected_tens_4[i], data_->tensor_4(q_point));
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
           	    // vector of local side quadrature points
        	    Range<EdgePoint> side_points = data_->side_eval->points(side, data_.get());
        	    for (EdgePoint side_p : side_points) {
                    EXPECT_EQ(expected_scalar[data_->computed_dh_cell_.elm_idx()], data_->scalar_field(side_p));
                    EXPECT_ARMA_EQ(expected_vector[i], data_->vector_field(side_p));
                    EXPECT_ARMA_EQ(expected_tens_1[i], data_->tensor_1(side_p));
                    EXPECT_ARMA_EQ(expected_tens_2[i], data_->tensor_2(side_p));
                    EXPECT_ARMA_EQ(expected_tens_3[i], data_->tensor_3(side_p));
                    EXPECT_ARMA_EQ(expected_tens_4[i], data_->tensor_4(side_p));
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
