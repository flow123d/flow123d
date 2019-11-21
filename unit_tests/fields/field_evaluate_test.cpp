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
        }

        // fields
        Field<3, FieldValue<3>::Scalar > scalar_field;
        Field<3, FieldValue<3>::VectorFixed > vector_field;
        Field<3, FieldValue<3>::TensorFixed > tensor_field;
    };

    FieldEval() {
        Profiler::initialize();
        data_ = std::make_shared<EqData>();
    }

    ~FieldEval() {
        std::cout << "-- destructor \n";
    }

    std::shared_ptr<EqData> data_;
};


/*TEST_F(FieldEval, eval_3d) {
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

    Profiler::initialize();
    PetscInitialize(0,PETSC_NULL,PETSC_NULL,PETSC_NULL);

    Mesh * mesh = mesh_full_constructor("{mesh_file=\"mesh/simplest_cube.msh\"}");
    std::shared_ptr<DOFHandlerMultiDim> dh = std::make_shared<DOFHandlerMultiDim>(*mesh);

    /// this can be done at initialization of the equation
	std::shared_ptr<EvalPoints> feval = std::make_shared<EvalPoints>();
    Quadrature *q_bulk = new QGauss(3, 2);
    Quadrature *q_side = new QGauss(2, 2);
    EvalSubset bulk_points = feval->add_bulk<3>(*q_bulk );
    EvalSubset side_points = feval->add_side<3>(*q_side );
    DHCellAccessor dh_cell(dh.get(), 3);

    std::cout << "Print bulk points:" << std::endl;
    for (auto p : bulk_points.points(dh_cell)) {
        std::cout << "--- bulk point:" << std::endl << p.loc_coords();
    }
    std::cout << "Print side points:" << std::endl;
    for (auto side_acc : dh_cell.side_range()) {
        std::cout << "- side idx: " << side_acc.side_idx() << ", permutation: " << side_acc.element()->permutation_idx( side_acc.side_idx() ) << std::endl;
        for ( auto p : side_points.points(side_acc) ) {
            std::cout << "--- side point" << std::endl << p.loc_coords();
        }
    }
  	std::cout << "----------- end \n";
}*/

TEST_F(FieldEval, evaluate) {
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

    Profiler::initialize();
    PetscInitialize(0,PETSC_NULL,PETSC_NULL,PETSC_NULL);

    Mesh * mesh = mesh_full_constructor("{mesh_file=\"mesh/simplest_cube.msh\"}");
    std::shared_ptr<DOFHandlerMultiDim> dh = std::make_shared<DOFHandlerMultiDim>(*mesh);

    // Asumme following types:
	std::shared_ptr<EvalPoints> feval = std::make_shared<EvalPoints>();
    Quadrature *q_bulk = new QGauss(3, 2);
    Quadrature *q_side = new QGauss(2, 2);
    EvalSubset mass_eval = feval->add_bulk<3>(*q_bulk );
    EvalSubset side_eval = feval->add_side<3>(*q_side );
    //EvalSubset this->ngh_side_eval;

    data_->cache_allocate(mass_eval);
    data_->cache_allocate(side_eval);

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
