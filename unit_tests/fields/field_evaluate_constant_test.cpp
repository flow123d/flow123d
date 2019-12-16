/*
 * field_evaluate_const_test.cpp
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
#include "fields/bc_field.hh"
#include "quadrature/quadrature.hh"
#include "quadrature/quadrature_lib.hh"
#include "fem/dofhandler.hh"
#include "fem/dh_cell_accessor.hh"
#include "mesh/mesh.h"
#include "mesh/sides.h"
#include "mesh/side_impl.hh"
#include "input/input_type.hh"
#include "input/accessors.hh"
#include "input/reader_to_storage.hh"
#include "system/sys_profiler.hh"


class FieldEvalConstantTest : public testing::Test {

public:
    class EqData : public FieldSet {
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

            for (unsigned int i=0; i<3; ++i)
                elm_cache_map_[i].init(i+1);
        }

        ElementCacheMap *get_element_cache_map(unsigned int dim) {
            return &elm_cache_map_[dim-1];
        }

        /// Add DHCellAccessor to appropriate ElementDataCache.
        void add_cell_to_cache(const DHCellAccessor &cell) {
        	elm_cache_map_[cell.dim()-1].add(cell);
        }


        // fields
        Field<3, FieldValue<3>::Scalar > scalar_field;
        Field<3, FieldValue<3>::VectorFixed > vector_field;
        Field<3, FieldValue<3>::TensorFixed > tensor_field;
        /// Element cache map of dimensions 1,2,3
        std::array< ElementCacheMap, 3 > elm_cache_map_;

    };

    FieldEvalConstantTest() {
        FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
        Profiler::initialize();
        PetscInitialize(0,PETSC_NULL,PETSC_NULL,PETSC_NULL);

        data_ = std::make_shared<EqData>();
        mesh_ = mesh_full_constructor("{mesh_file=\"mesh/cube_2x1.msh\"}");
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
    }

    std::shared_ptr<EqData> data_;
    Mesh * mesh_;
    std::shared_ptr<DOFHandlerMultiDim> dh_;
};


TEST_F(FieldEvalConstantTest, evaluate) {
    string eq_data_input = R"YAML(
    data:
      - region: 3D left
        time: 0.0
        scalar_field: !FieldConstant
          value: 0.5
        vector_field: [1, 2, 3]
        tensor_field: [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
      - region: 3D right
        time: 0.0
        scalar_field: !FieldConstant
          value: 1.5
        vector_field: [4, 5, 6]
        tensor_field: [2.1, 2.2, 2.3, 2.4, 2.5, 2.6]
    )YAML";
	this->read_input(eq_data_input);

    // Asumme following types:
	std::shared_ptr<EvalPoints> feval = std::make_shared<EvalPoints>();
    Quadrature *q_bulk = new QGauss(3, 2);
    Quadrature *q_side = new QGauss(2, 2);
    std::shared_ptr<EvalSubset> mass_eval = feval->add_bulk<3>(*q_bulk );
    std::shared_ptr<EvalSubset> side_eval = feval->add_side<3>(*q_side );
    //std::shared_ptr<EvalSubset> this->ngh_side_eval;

    data_->cache_allocate(mass_eval);
    data_->cache_allocate(side_eval);

    std::vector<unsigned int> cell_idx = {3, 4, 5, 10};
    std::vector<double>       expected_scalar = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5};
    std::vector<arma::vec3>   expected_vector = {{1, 2, 3}, {1, 2, 3}, {1, 2, 3}, {4, 5, 6}};
    std::vector<arma::mat33>  expected_tensor = {{0.1, 0.2, 0.3, 0.2, 0.4, 0.5, 0.3, 0.5, 0.6}, {0.1, 0.2, 0.3, 0.2, 0.4, 0.5, 0.3, 0.5, 0.6},
                                                 {0.1, 0.2, 0.3, 0.2, 0.4, 0.5, 0.3, 0.5, 0.6}, {2.1, 2.2, 2.3, 2.2, 2.4, 2.5, 2.3, 2.5, 2.6}};
    for (unsigned int i=0; i<cell_idx.size(); ++i) {
        DHCellAccessor dh_cell(dh_.get(), cell_idx[i]);  // element ids stored to cache: (3 -> 2,3,4), (4 -> 3,4,5,10), (5 -> 0,4,5,11), (10 -> 4,9,10,11)
        data_->add_cell_to_cache(dh_cell);
        for (DHCellSide side : dh_cell.side_range()) {
            for(DHCellSide el_ngh_side : side.edge_sides()) {
                data_->add_cell_to_cache( el_ngh_side.cell() );
    	    }
        }
        data_->cache_update(*data_->get_element_cache_map(3), mesh_);

        DHCellAccessor cache_cell = this->data_->elm_cache_map_[dh_cell.dim()-1](dh_cell);

        // Bulk integral, no sides, no permutations.
        for(BulkPoint q_point: mass_eval->points(cache_cell)) {
            EXPECT_EQ(expected_scalar[cache_cell.elm_idx()], data_->scalar_field(q_point)[0]);
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
        for (DHCellSide side : cache_cell.side_range()) {
        	for(DHCellSide el_ngh_side : side.edge_sides()) {
        		el_ngh_side.cell() = this->data_->elm_cache_map_[el_ngh_side.cell().dim()-1](el_ngh_side.cell());
           	    // vector of local side quadrature points in the correct side permutation
        	    Range<SidePoint> side_points = side_eval->points(side);
        	    for (SidePoint side_p : side_points) {
                    EXPECT_EQ(expected_scalar[cache_cell.elm_idx()], data_->scalar_field(side_p)[0]);
                    EXPECT_ARMA_EQ(expected_vector[i], data_->vector_field(side_p));
                    EXPECT_ARMA_EQ(expected_tensor[i], data_->tensor_field(side_p));
                    SidePoint ngh_p = side_p.permute(el_ngh_side);
                    EXPECT_EQ(expected_scalar[el_ngh_side.cell().elm_idx()], data_->scalar_field(ngh_p)[0]);
        	        //loc_mat += cross_section(side_p) * sigma(side_p) *
        		    //    (conc.base_value(side_p) * velocity(side_p)
        		    //    + conc.base_value(ngh_p) * velocity(ngh_p)) * side_p.normal() / 2;
                }
            }
        }
    }

}
