/*
 * field_evaluate_formula_test.cpp
 *
 *  Created on: Dec 03, 2019
 *      Author: David Flanderka
 *
 *  Tests evaluation of FieldFormula
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


class FieldEvalFormulaTest : public testing::Test {

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
            this->cache_reallocate(*this);
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

        void reallocate_cache() {
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

    FieldEvalFormulaTest() {
        FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
        Profiler::instance();
        PetscInitialize(0,PETSC_NULL,PETSC_NULL,PETSC_NULL);

        data_ = std::make_shared<EqData>();
        data_->add_coords_field();
        mesh_ = mesh_full_constructor("{ mesh_file=\"mesh/cube_2x1.msh\", optimize_mesh=false }");
        dh_ = std::make_shared<DOFHandlerMultiDim>(*mesh_);
    }

    ~FieldEvalFormulaTest() {}

    static Input::Type::Record & get_input_type() {
        return IT::Record("SomeEquation","")
                .declare_key("data", IT::Array(
                        IT::Record("SomeEquation_Data", FieldCommon::field_descriptor_record_description("SomeEquation_Data") )
                        .copy_keys( FieldEvalFormulaTest::EqData().make_field_descriptor_type("SomeEquation") )
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
        data_->set_dependency();
        data_->reallocate_cache();
    }


    std::shared_ptr<EqData> data_;
    Mesh * mesh_;
    std::shared_ptr<DOFHandlerMultiDim> dh_;
};


TEST_F(FieldEvalFormulaTest, evaluate) {
    string eq_data_input = R"YAML(
    data:
      - region: 3D left
        time: 0.0
        scalar_field: !FieldFormula
          value: x
        vector_field: !FieldFormula
          value: [x, 2*x, 0.5]
        tensor_field: !FieldFormula
          value: [x, 0.2, 0.3, 0.4, 0.5, 0.6]
      - region: 3D right
        time: 0.0
        scalar_field: !FieldFormula
          value: y
        vector_field:  !FieldFormula
          value: [y, 2*y, 0.5]
        tensor_field: !FieldFormula
          value: [y, 2.2, 2.3, 2.4, 2.5, 2.6]
    )YAML";
	this->read_input(eq_data_input);

    std::vector<unsigned int> cell_idx = {3, 4, 5, 9};
    std::vector<double> dbl_scalar = { -1.447213595499958094, -0.552786404500041906, -0.276393202250021008, -1.170820393249937030, -0.723606797749979047,
                                       +0.170820393249937058, -1.666666666666666741, -1.333333333333333481, -0.333333333333333315, -0.666666666666666519,
                                       +0.000000000000000000, +0.333333333333333481, -1.000000000000000000, +0.666666666666666741 };
    std::vector< std::vector<uint> > expected_scalar = { {0,1,0,1,9,9,6,6,6,6,9,9,6,6,6,6,7,7,8,8,8,8,6,6,9,9,6,6,8,8,8,8,7,7,8,8,8,8,7,7},
            {2,2,3,2,8,8,8,8,7,7,8,8,8,8,7,7,10,10,10,10,10,10,10,9,10,9,10,11,8,8,8,8,7,7,8,8,8,8,7,7,8,8,8,8,7,7},
			{2,2,3,2,8,8,8,8,7,7,8,8,8,8,7,7,10,10,10,10,10,10,10,8,10,13,10,13,8,8,8,8,7,7,8,8,7,7,8,8,8,8,7,7,8,8},
			{4,4,4,5,9,9,9,9,11,11,9,9,9,9,11,11,9,9,9,9,11,11,9,9,9,9,11,11,12,12,12,12,12,12,11,11,9,9,9,9} };
    std::vector<arma::mat33>  expected_tensor = {{0.1, 0.2, 0.3, 0.2, 0.4, 0.5, 0.3, 0.5, 0.6}, {0.1, 0.2, 0.3, 0.2, 0.4, 0.5, 0.3, 0.5, 0.6},
                                                 {0.1, 0.2, 0.3, 0.2, 0.4, 0.5, 0.3, 0.5, 0.6}, {2.1, 2.2, 2.3, 2.2, 2.4, 2.5, 2.3, 2.5, 2.6}};
    for (uint i=0; i<cell_idx.size(); ++i) {
        uint test_point = 0; // index to expected vals
    	data_->start_elements_update();
    	data_->computed_dh_cell_ = DHCellAccessor(dh_.get(), cell_idx[i]);  // element ids stored to cache: (3 -> 2,3,4), (4 -> 3,4,5,10), (5 -> 0,4,5,11), (10 -> 8,9,10)
        data_->update_cache();

        // Bulk integral, no sides, no permutations.
        for( BulkPoint q_point: data_->mass_eval->points(data_->computed_dh_cell_, data_.get()) ) {
            EXPECT_DOUBLE_EQ( dbl_scalar[ expected_scalar[i][test_point] ], data_->scalar_field(q_point));
            arma::vec3 expected_vector;
            expected_vector(0) = dbl_scalar[ expected_scalar[i][test_point] ];
            expected_vector(1) = 2*expected_vector(0);
            expected_vector(2) = 0.5;
            EXPECT_ARMA_EQ(expected_vector, data_->vector_field(q_point));
            auto exp_tensor = expected_tensor[i];
            exp_tensor(0,0) = dbl_scalar[ expected_scalar[i][test_point] ];
            EXPECT_ARMA_EQ(exp_tensor, data_->tensor_field(q_point));
            test_point++;
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
                    EXPECT_DOUBLE_EQ( dbl_scalar[ expected_scalar[i][test_point] ], data_->scalar_field(side_p));
                    arma::vec3 expected_vector;
                    expected_vector(0) = dbl_scalar[ expected_scalar[i][test_point] ];
                    expected_vector(1) = 2*expected_vector(0);
                    expected_vector(2) = 0.5;
                    EXPECT_ARMA_EQ(expected_vector, data_->vector_field(side_p));
                    auto exp_tensor = expected_tensor[i];
                    exp_tensor(0,0) = dbl_scalar[ expected_scalar[i][test_point] ];
                    EXPECT_ARMA_EQ(exp_tensor, data_->tensor_field(side_p));
                    test_point++;
                    EdgePoint ngh_p = side_p.point_on(el_ngh_side);
                    EXPECT_DOUBLE_EQ( dbl_scalar[ expected_scalar[i][test_point] ], data_->scalar_field(ngh_p));
                    test_point++;
        	        //loc_mat += cross_section(side_p) * sigma(side_p) *
        		    //    (conc.base_value(side_p) * velocity(side_p)
        		    //    + conc.base_value(ngh_p) * velocity(ngh_p)) * side_p.normal() / 2;
                }
            }
        }
    }

}
