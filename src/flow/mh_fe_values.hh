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

#include "fields/field_base.hh"
#include "fields/field_values.hh"

typedef Field<3, FieldValue<3>::TensorFixed > FieldType;
typedef Field<3, FieldValue<3>::Scalar > FieldType_Scalar;

//declared in system/math_fce.h
//typedef double SmallVec2_t[2];
//typedef SmallVec2_t *SmallMtx2;

/**
 * Temporary class to remove MH calculations (basis functions and leading local matrix) from geometrical mesh.
 */
class MHFEValues {
public:
    MHFEValues();
    ~MHFEValues();
    void update(ElementFullIter ele, FieldType &cond_anisothropy, FieldType_Scalar &cross_section);
    double * local_matrix();
    double * inv_local_matrix();

    /**
     * Temporary hack: returns value of shape function on element 'ele' and its 'face' in 'point' given in global coordinate system.
     */
    arma::vec3 RT0_value(ElementFullIter ele, arma::vec3 point, unsigned int face);

private:
    void local_matrix_line(ElementFullIter ele, FieldType &cond_anisothropy, FieldType_Scalar &cross_section);
    void local_matrix_triangle(ElementFullIter ele, FieldType &cond_anisothropy, FieldType_Scalar &cross_section);
    void local_matrix_tetrahedron(ElementFullIter ele, FieldType &cond_anisothropy, FieldType_Scalar &cross_section);


    void node_coordinates_triangle( ElementFullIter ele, double nod[ 3 ][ 2 ] );
    void side_midpoint_triangle( double nod[ 3 ][ 2 ], double midpoint[ 3 ][ 2 ] );
    void basis_functions_triangle( double nod[ 3 ][ 2 ], double alfa[],double beta[], double gama[] );
    void bas_func_0_triangle( double x0, double y0,
                              double x1, double y1,
                              double x2, double y2,
                              double *alfa, double *beta, double *gama);
    void calc_polynom_triangle( double al_i, double be_i, double al_j, double be_j, arma::mat::fixed<2,2> a, double poly[] );
    double polynom_value_triangle( double poly[], double point[] );
    void basis_functions_tetrahedron( ElementFullIter ele, double alfa[], double beta[],double gama[], double delta[], FieldType_Scalar &cross_section);
    void calc_polynom_tetrahedron( double al_i, double be_i, double ga_i,
                                   double al_j, double be_j, double ga_j,
                                   arma::mat::fixed<3,3> a, double poly[] );
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
