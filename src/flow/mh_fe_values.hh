/*
 * mh_fe_values.hh
 *
 *  Created on: Jul 2, 2012
 *      Author: jb
 */

#ifndef MH_FE_VALUES_HH_
#define MH_FE_VALUES_HH_

#include "mesh/mesh_types.hh"
#include <armadillo>


/**
 * Temporary class to remove MH calculations (basis functions and leading local matrix) from geometrical mesh.
 */
class MHFEValues {
public:
    MHFEValues();
    ~MHFEValues();
    void update(ElementFullIter ele);
    double * local_matrix();
    double * inv_local_matrix();

    arma::vec3 RT0_value(ElementFullIter ele, arma::vec3 point, unsigned int face);

private:
    void local_matrix_line(ElementFullIter ele );
    void local_matrix_triangle(ElementFullIter ele );
    void local_matrix_tetrahedron(ElementFullIter ele );


    void node_coordinates_triangle( ElementFullIter ele, double nod[ 3 ][ 2 ] );
    void side_midpoint_triangle( double nod[ 3 ][ 2 ], double midpoint[ 3 ][ 2 ] );
    void basis_functions_triangle( double nod[ 3 ][ 2 ], double alfa[],double beta[], double gama[] );
    void bas_func_0_triangle( double x0, double y0,
                              double x1, double y1,
                              double x2, double y2,
                              double *alfa, double *beta, double *gama);
    void calc_polynom_triangle( double al_i, double be_i, double al_j, double be_j,SmallMtx2 a, double poly[] );
    double polynom_value_triangle( double poly[], double point[] );
    void basis_functions_tetrahedron( ElementFullIter ele, double alfa[], double beta[],double gama[], double delta[]);
    void calc_polynom_tetrahedron( double al_i, double be_i, double ga_i,
                                   double al_j, double be_j, double ga_j,
                                   double a[ 3 ][ 3 ], double poly[] );
    double polynom_integral_tetrahedron( ElementFullIter ele, double poly[] );

    double * loc_matrix_;
    double * inv_loc_matrix_;


    // Parameters of the basis functions
    double   *bas_alfa;      // Parameters alfa
    double   *bas_beta;      // Parameters beta
    double   *bas_gama;      // Parameters gama
    double   *bas_delta;      // Parameters delta
};


#endif /* MH_FE_VALUES_HH_ */
