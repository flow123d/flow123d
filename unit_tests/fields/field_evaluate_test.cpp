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
#include "system/sys_profiler.hh"


class FieldEval : public testing::Test {

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

    FieldEval() {
        Profiler::initialize();
        data_ = std::make_shared<EqData>();
    }

    ~FieldEval() {}

    std::shared_ptr<EqData> data_;
};


TEST_F(FieldEval, eval_3d) {
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

    Profiler::initialize();
    PetscInitialize(0,PETSC_NULL,PETSC_NULL,PETSC_NULL);

    Mesh * mesh = mesh_full_constructor("{mesh_file=\"mesh/simplest_cube.msh\"}");
    std::shared_ptr<DOFHandlerMultiDim> dh = std::make_shared<DOFHandlerMultiDim>(*mesh);

    /// this can be done at initialization of the equation
	std::shared_ptr<EvalPoints> feval = std::make_shared<EvalPoints>();
    Quadrature *q_bulk = new QGauss(3, 2);
    Quadrature *q_side = new QGauss(2, 2);
    std::shared_ptr<EvalSubset> bulk_points = feval->add_bulk<3>(*q_bulk );
    std::shared_ptr<EvalSubset> side_points = feval->add_side<3>(*q_side );
    DHCellAccessor dh_cell(dh.get(), 3);

    {
        // Test of bulk local points
    	std::vector<arma::vec3> expected_vals = {{0.138196601125010504, 0.138196601125010504, 0.138196601125010504},
    			                                 {0.138196601125010504, 0.138196601125010504, 0.585410196624968515},
												 {0.138196601125010504, 0.585410196624968515, 0.138196601125010504},
												 {0.585410196624968515, 0.138196601125010504, 0.138196601125010504}};
    	unsigned int i=0; // iter trought expected_vals
    	for (auto p : bulk_points->points(dh_cell)) {
            EXPECT_ARMA_EQ(p.loc_coords(), expected_vals[i]);
			++i;
        }
    }
    {
        // Test of side local points
        std::vector< std::vector<arma::vec3> > expected_vals(4);
        expected_vals[0] = { {0.166666666666666657, 0.166666666666666657, 0.0},
                             {0.666666666666666741, 0.166666666666666657, 0.0},
                             {0.166666666666666657, 0.666666666666666741, 0.0} };
        expected_vals[1] = { {0.166666666666666657, 0.0, 0.166666666666666657},
                             {0.166666666666666657, 0.0, 0.666666666666666741},
                             {0.666666666666666741, 0.0, 0.166666666666666657} };
        expected_vals[2] = { {0.0, 0.166666666666666657, 0.166666666666666657},
	                         {0.0, 0.166666666666666657, 0.666666666666666741},
                             {0.0, 0.666666666666666741, 0.166666666666666657} };
        expected_vals[3] = { {0.666666666666666741, 0.166666666666666657, 0.166666666666666657},
                             {0.166666666666666657, 0.166666666666666657, 0.666666666666666741},
                             {0.166666666666666657, 0.666666666666666741, 0.166666666666666657} };
        unsigned int i_side=0, i_point; // iter trought expected_vals
        for (auto side_acc : dh_cell.side_range()) {
        	i_point=0;
            for ( auto p : side_points->points(side_acc) ) {
            	EXPECT_ARMA_EQ(p.loc_coords(), expected_vals[i_side][i_point]);
                ++i_point;
            }
            ++i_side;
        }
    }
}

/*
 * Prepared test for further development
 */
TEST_F(FieldEval, evaluate) {
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

    Profiler::initialize();
    PetscInitialize(0,PETSC_NULL,PETSC_NULL,PETSC_NULL);

    Mesh * mesh = mesh_full_constructor("{mesh_file=\"mesh/simplest_cube.msh\"}");
    data_->set_mesh(*mesh);
    std::shared_ptr<DOFHandlerMultiDim> dh = std::make_shared<DOFHandlerMultiDim>(*mesh);

    // Asumme following types:
	std::shared_ptr<EvalPoints> feval = std::make_shared<EvalPoints>();
    Quadrature *q_bulk = new QGauss(3, 2);
    Quadrature *q_side = new QGauss(2, 2);
    std::shared_ptr<EvalSubset> mass_eval = feval->add_bulk<3>(*q_bulk );
    std::shared_ptr<EvalSubset> side_eval = feval->add_side<3>(*q_side );
    //std::shared_ptr<EvalSubset> this->ngh_side_eval;

    data_->cache_allocate(mass_eval, data_->get_element_cache_map(3));
    data_->cache_allocate(side_eval, data_->get_element_cache_map(3));

    //DHCellAccessor cache_cell = this->element_cache_map(cell);
    DHCellAccessor cache_cell(dh.get(), 4);  // element ids store to cache: (3 -> 3,4), (4 -> 3,4,5), (5 -> 4,5)
    data_->add_cell_to_cache(cache_cell);
    for (DHCellSide side : cache_cell.side_range()) {
    	for(DHCellSide el_ngh_side : side.edge_sides()) {
    	    data_->add_cell_to_cache( el_ngh_side.cell() );
    	}
    }
    //data_->cache_update(data_->get_element_cache_map(3));

    //...
    /*DHCellAccessor cache_cell = this->element_cache_map(cell);
    // Bulk integral, no sides, no permutations.
    for(BulkPoint q_point: this->mass_eval.points(cache_cell)) {
        // Extracting the cached values.
        double cs = cross_section(q_point);

        // Following would be nice to have. Not clear how to
        // deal with more then single element as fe_values have its own cache that has to be updated.
        auto base_fn_grad = presssure_field_fe.base_value(q_point);
    loc_matrix += outer_product((cs * base_fn_grad),  base_fn_grad)
    } */

    // Side integrals.
    // FieldFE<..> conc;
    /*for (DHCellSide side : cache_cell.side_range()) {
    	for(DHCellSide el_ngh_side : side.edge_sides()) {
       	    // vector of local side quadrature points in the correct side permutation
    	    Range<SidePoint> side_points = this->side_eval.points(side)
    	    for (SidePoint p : side_points) {
    	    	ngh_p = p.permute(el_ngh_side);
    	        loc_mat += cross_section(p) * sigma(p) *
    		    (conc.base_value(p) * velocity(p)
    		    + conc.base_value(ngh_p) * velocity(ngh_p)) * p.normal() / 2;
            }
        }
    }*/
    //std::cout << "----------- end \n";

}
