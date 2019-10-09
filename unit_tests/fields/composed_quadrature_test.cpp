#define TEST_USE_MPI
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>
#include <mesh_constructor.hh>

#include "fields/eval_points.hh"
#include "fields/eval_subset.hh"
#include "fields/field_value_cache.hh"
#include "fields/field_values.hh"
#include "quadrature/quadrature.hh"
#include "quadrature/quadrature_lib.hh"
#include "fem/dofhandler.hh"
#include "fem/dh_cell_accessor.hh"
#include "mesh/mesh.h"
#include "mesh/sides.h"
#include "mesh/side_impl.hh"
#include "system/sys_profiler.hh"


TEST(ComposedQuadratureTest, eval_3d) {
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

    Profiler::initialize();
    PetscInitialize(0,PETSC_NULL,PETSC_NULL,PETSC_NULL);

    Mesh * mesh = mesh_full_constructor("{mesh_file=\"mesh/simplest_cube.msh\"}");
    std::shared_ptr<DOFHandlerMultiDim> dh = std::make_shared<DOFHandlerMultiDim>(*mesh);

    /// this can be done at initialization of the equation
	EvalPoints feval;
    Quadrature<3> *q_bulk = new QGauss<3>(2);
    Quadrature<2> *q_side = new QGauss<2>(2);
    EvalSubset bulk_points = feval.add_bulk<3>(*q_bulk );
    EvalSubset side_points = feval.add_side<3>(*q_side );
    DHCellAccessor dh_cell(dh.get(), 3);
    /// this is part of assembly process
    //for (auto cell : dh->own_range()) {
    //    feval.reinit(cell.elm());
        std::cout << "Print bulk points:" << std::endl;
        for (auto p : bulk_points.points(dh_cell)) {
        	std::cout << "--- bulk point:" << std::endl << p.loc_coords();
            //double bulk_expr = cross_section.get_value(p) * conductivity.get_value(p);
        }
        std::cout << "Print side points:" << std::endl;
        for (auto side_acc : dh_cell.side_range()) {
            std::cout << "- side idx: " << side_acc.side_idx() << ", permutation: " << side_acc.element()->permutation_idx( side_acc.side_idx() ) << std::endl;
            for ( auto p : side_points.points(side_acc) ) {
                std::cout << "--- side point" << std::endl << p.loc_coords();
                //double side_expr = cross_section.get_value(p) * sigma.get_value(p);
            }
        }
    //}
  	std::cout << "----------- end \n";

  	FieldValueCache<FieldValue<0>::Scalar> scalar_cache(side_points.eval_points());
  	FieldValueCache<FieldValue<3>::VectorFixed> vector_cache(side_points.eval_points());
  	FieldValueCache<FieldValue<3>::TensorFixed> tensor_cache(side_points.eval_points());
  	scalar_cache.mark_used(side_points);
  	vector_cache.mark_used(side_points);
  	tensor_cache.mark_used(side_points);
}
