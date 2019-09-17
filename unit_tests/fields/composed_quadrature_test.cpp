#define TEST_USE_MPI
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>

#include "fields/composed_quadrature.hh"
#include "fields/point_sets.hh"
#include "quadrature/quadrature.hh"
#include "quadrature/quadrature_lib.hh"
#include "mesh/sides.h"


TEST(ComposedQuadratureTest, eval_3d) {
    /// this can be done at initialization of the equation
	ComposedQuadrature<3> feval;
    Quadrature<3> *q_bulk = new QGauss<3>(2);
    Quadrature<2> *q_side = new QGauss<2>(2);
    BulkSubQuad<3> bulk_points = feval.add_bulk_quad(*q_bulk );
    SideSubQuad<3> side_points = feval.add_side_quad(*q_side );
    /// this is part of assembly process
    //for (auto cell : dh->own_range()) {
    //    feval.reinit(cell.elm());
        for (auto p : bulk_points.points()) {
        	std::cout << "--- bulk point:" << std::endl << p.loc_coords();
            //double bulk_expr = cross_section.get_value(p) * conductivity.get_value(p);
        }
        for (Side side; side.side_idx()<4; side.inc()) {
            for ( auto p : side_points.points(side) ) {
            	std::cout << "--- side point " << side.side_idx() << std::endl << p.loc_coords();
                //double side_expr = cross_section.get_value(p) * sigma.get_value(p);
            }
        }
    //}
  	std::cout << "----------- end \n";
}
