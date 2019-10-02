#define TEST_USE_MPI
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>

#include "fields/eval_points.hh"
#include "fields/eval_subset.hh"
#include "fields/field_value_cache.hh"
#include "fields/field_values.hh"
#include "quadrature/quadrature.hh"
#include "quadrature/quadrature_lib.hh"
#include "mesh/sides.h"
#include "mesh/side_impl.hh"


TEST(ComposedQuadratureTest, eval_3d) {
    /// this can be done at initialization of the equation
	EvalPoints feval;
    Quadrature<3> *q_bulk = new QGauss<3>(2);
    Quadrature<2> *q_side = new QGauss<2>(2);
    EvalSubset bulk_points = feval.add_bulk<3>(*q_bulk );
    EvalSubset side_points = feval.add_side<3>(*q_side );
    bulk_points.print_bulk_points();
    side_points.print_side_points(0);
  	std::cout << "----------- end \n";
  	FieldValueCache<FieldValue<0>::Scalar> scalar_cache(feval);
  	FieldValueCache<FieldValue<3>::VectorFixed> vector_cache(feval);
  	FieldValueCache<FieldValue<3>::TensorFixed> tensor_cache(feval);
}
