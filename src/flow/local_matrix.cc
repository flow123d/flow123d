/*!
 *
 * Copyright (C) 2007 Technical University of Liberec.  All rights reserved.
 *
 * Please make a following refer to Flow123d on your project site if you use the program for any purpose,
 * especially for academic research:
 * Flow123d, Research Centre: Advanced Remedial Technologies, Technical University of Liberec, Czech Republic
 *
 * This program is free software; you can redistribute it and/or modify it under the terms
 * of the GNU General Public License version 3 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 021110-1307, USA.
 *
 *
 * $Id$
 * $Revision$
 * $LastChangedBy$
 * $LastChangedDate$
 *
 * @file
 * @brief  Compute local matrix
 *
 */

#include "system/system.hh"
#include "system/math_fce.h"
#include "mesh/mesh.h"
#include "flow/local_matrix.h"
#include "materials.hh"

static void local_matrix(ElementIter  ele);
static void local_matrix_triangle(ElementIter  ele);
static void node_coordinates_triangle(ElementIter  ele,double[3][2]);
static void side_midpoint_triangle(double[3][2],double[3][2]);
static void basis_functions_triangle(double[3][2],double[],double[],double[]);
static void bas_func_0_triangle(double,double,double,double,double,double,double*,double*,double*);
static void calc_polynom_triangle(double,double,double,double,double[2][2],double[]);
static double polynom_value_triangle(double[],double[]);
static void local_matrix_line(ElementIter );
static void local_matrix_tetrahedron(ElementIter );
static void basis_functions_tetrahedron(ElementIter ,double[],double[],double[],double[]);
static void calc_polynom_tetrahedron(double,double,double,double,double,double,double[3][3],double[]);
static double polynom_integral_tetrahedron(ElementIter ,double[]);

//=============================================================================
// CALCULATE LOCAL MATRICES OF ALL ELEMENTS OF THE MESH
//=============================================================================
void local_matrices_mh(Mesh* mesh)
{
	ElementIter ele;

	xprintf( Msg, "Calculating local matrices... ")/*orig verb 2*/;
	ASSERT(NONULL( mesh ),"No mesh for problem\n");

	FOR_ELEMENTS(mesh, ele )
		local_matrix( ele );
}
//=============================================================================
// CALCULATE LOCAL MARIX FOR ONE ELEMENT
//=============================================================================
void local_matrix( ElementIter ele )
{
	ASSERT(!( ele == NULL ),"NULL as argument of function local_matrix()\n");
	ele->loc = (double *) xmalloc( ele->n_sides * ele->n_sides*sizeof(double));
	switch( ele->type ) {
		case LINE:
			local_matrix_line( ele );
			break;
		case TRIANGLE:
			local_matrix_triangle( ele );
			break;
		case TETRAHEDRON:
			local_matrix_tetrahedron( ele );
			break;
	}
}
//=============================================================================
// CALCULATE LOCAL MATRIX FOR TRIANGULAR ELEMENT
//=============================================================================
void local_matrix_triangle( ElementIter ele )
{
	double midpoint[ 3 ][ 2 ]; // Midpoints of element's sides
	double alfa[ 3 ];	//
	double beta[ 3 ];       //  |- Parametrs of basis functions
	double gama[ 3 ];       // /
	double poly[ 6 ];       // Koeficients of polynomial
	double p[ 3 ];         	// Polynomial's values in midpoints
	int i, j;               // Loops' counters
	double nod_coor[ 3 ][ 2 ]; // Coordinates of nodes;
	SmallMtx3 loc=(SmallMtx3)(ele->loc);

	node_coordinates_triangle( ele, nod_coor );
   	side_midpoint_triangle( nod_coor, midpoint );
	basis_functions_triangle( nod_coor, alfa, beta, gama );
	for( i = 0; i < 3; i++ )
		for( j = i; j < 3; j++ ) {
			calc_polynom_triangle( alfa[ i ], beta[ i ],
			          	       alfa[ j ], beta[ j ],
			          	     (SmallMtx2)(ele->a), poly );
			p[ 0 ] = polynom_value_triangle( poly, midpoint[ 0 ] );
			p[ 1 ] = polynom_value_triangle( poly, midpoint[ 1 ] );
			p[ 2 ] = polynom_value_triangle( poly, midpoint[ 2 ] );
			loc[i][j] = gama[i] * gama[j] * ele->measure() *
					( p[0] + p[1] + p[2] ) / (3.0 * ele->material->size);
			loc[j][i] = loc[i][j];
		}
        for(i = 0; i < 3; i++) {
                ele->bas_alfa[ i ] = alfa[ i ];
                ele->bas_beta[ i ] = beta[ i ];
                ele->bas_gama[ i ] = gama[ i ];
        }
	//check_local( ele );
}
//=============================================================================
//
//=============================================================================
void node_coordinates_triangle( ElementIter ele, double nod[ 3 ][ 2 ] )
{
	double u[ 3 ], v[ 3 ];	// vectors of sides 0 and 2
	double t[ 3 ];		// u normalized (tangenta)
	double n[ 3 ];		// normal
	double b[ 3 ];		// binormal

	u[ 0 ] = ele->node[ 1 ]->getX() - ele->node[ 0 ]->getX();
	u[ 1 ] = ele->node[ 1 ]->getY() - ele->node[ 0 ]->getY();
	u[ 2 ] = ele->node[ 1 ]->getZ() - ele->node[ 0 ]->getZ();
	v[ 0 ] = ele->node[ 2 ]->getX() - ele->node[ 0 ]->getX();
	v[ 1 ] = ele->node[ 2 ]->getY() - ele->node[ 0 ]->getY();
	v[ 2 ] = ele->node[ 2 ]->getZ() - ele->node[ 0 ]->getZ();

	vector_product( u, v, n );
	t[ 0 ] = u[ 0 ];
	t[ 1 ] = u[ 1 ];
	t[ 2 ] = u[ 2 ];
	normalize_vector( t );
	normalize_vector( n );
	vector_product( n, u, b );
	normalize_vector( b );

	nod[ 0 ][ 0 ] = 0.0;
	nod[ 0 ][ 1 ] = 0.0;
	nod[ 1 ][ 0 ] = scalar_product( t, u );
	nod[ 1 ][ 1 ] = scalar_product( b, u );
	nod[ 2 ][ 0 ] = scalar_product( t, v );
	nod[ 2 ][ 1 ] = scalar_product( b, v );
}
//=============================================================================
// CALCULATE LOCAL MATRIX FOR TRIANGLE ELEMENT
//=============================================================================
void side_midpoint_triangle( double nod[ 3 ][ 2 ], double midpoint[ 3 ][ 2 ] )
{
	midpoint[ 0 ][ 0 ] = ( nod[ 0 ][ 0 ] + nod[ 1 ][ 0 ] ) / 2.0;
	midpoint[ 0 ][ 1 ] = ( nod[ 0 ][ 1 ] + nod[ 1 ][ 1 ] ) / 2.0;
	midpoint[ 1 ][ 0 ] = ( nod[ 1 ][ 0 ] + nod[ 2 ][ 0 ] ) / 2.0;
	midpoint[ 1 ][ 1 ] = ( nod[ 1 ][ 1 ] + nod[ 2 ][ 1 ] ) / 2.0;
	midpoint[ 2 ][ 0 ] = ( nod[ 2 ][ 0 ] + nod[ 0 ][ 0 ] ) / 2.0;
	midpoint[ 2 ][ 1 ] = ( nod[ 2 ][ 1 ] + nod[ 0 ][ 1 ] ) / 2.0;
}
//=============================================================================
// GENERATE BASIS FUNCTIONS
//=============================================================================
void basis_functions_triangle( double nod[ 3 ][ 2 ], double alfa[],
		      	       double beta[], double gama[] )
{
        bas_func_0_triangle( nod[ 0 ][ 0 ], nod[ 0 ][ 1 ],
	                     nod[ 1 ][ 0 ], nod[ 1 ][ 1 ],
         	             nod[ 2 ][ 0 ], nod[ 2 ][ 1 ],
                  	     alfa, beta, gama );
        bas_func_0_triangle( nod[ 1 ][ 0 ], nod[ 1 ][ 1 ],
			     nod[ 2 ][ 0 ], nod[ 2 ][ 1 ],
	                     nod[ 0 ][ 0 ], nod[ 0 ][ 1 ],
         	             alfa + 1, beta + 1, gama + 1 );
        bas_func_0_triangle( nod[ 2 ][ 0 ], nod[ 2 ][ 1 ],
                  	     nod[ 0 ][ 0 ], nod[ 0 ][ 1 ],
	                     nod[ 1 ][ 0 ], nod[ 1 ][ 1 ],
         	             alfa + 2, beta + 2, gama + 2 );
}
//=============================================================================
// GENERATE BASIS FUNCTION FOR SIDE 0
//=============================================================================
void bas_func_0_triangle( double x0, double y0,
                          double x1, double y1,
                          double x2, double y2,
                          double *alfa, double *beta, double *gama)
{
	*alfa = x2;
	*beta = y2;
	*gama = 1.0 / fabs( ( x0  - x2 ) * ( y1 - y0 ) + ( y0 - y2 ) * ( x0 - x1 ) );
}
//=============================================================================
// CALCULATE POLYNOM OF SCALAR PRODUCT
//=============================================================================
void calc_polynom_triangle( double al_i, double be_i, double al_j, double be_j,
	                    SmallMtx2 a, double poly[] )
{
        poly[ 0 ] =   a[ 0 ][ 0 ] * al_i * al_j +
		      a[ 0 ][ 1 ] * be_i * al_j +
		      a[ 1 ][ 0 ] * al_i * be_j +
                      a[ 1 ][ 1 ] * be_i * be_j;
        poly[ 1 ] = ( a[ 0 ][ 0 ] * al_i +
		      a[ 0 ][ 0 ] * al_j +
		      a[ 0 ][ 1 ] * be_i +
		      a[ 1 ][ 0 ] * be_j ) * -1.0;
        poly[ 2 ] = ( a[ 1 ][ 1 ] * be_i +
		      a[ 1 ][ 1 ] * be_j +
		      a[ 1 ][ 0 ] * al_i +
		      a[ 0 ][ 1 ] * al_j ) * -1.0;
        poly[ 3 ] =   a[ 0 ][ 0 ];
        poly[ 4 ] =   a[ 0 ][ 1 ] + a[ 1 ][ 0 ];
        poly[ 5 ] =   a[ 1 ][ 1 ];
}
//=============================================================================
// CALCULATE VALUE OF POLYNOM IN GIVEN POINT
//=============================================================================
double polynom_value_triangle( double poly[], double point[] )
{
	double rc;

	rc = poly[ 0 ] +
	     poly[ 1 ] * point[ 0 ] +
	     poly[ 2 ] * point[ 1 ] +
	     poly[ 3 ] * point[ 0 ] * point[ 0 ] +
	     poly[ 4 ] * point[ 0 ] * point[ 1 ] +
	     poly[ 5 ] * point[ 1 ] * point[ 1 ];
	return rc;
}
//=============================================================================
// CALCULATE LOCAL MARIX FOR LINEAR ELEMENT
//=============================================================================
void local_matrix_line( ElementIter ele )
{
	double 		val;
	SmallMtx2 loc=(SmallMtx2)(ele->loc);

	val= (* ele->a) * ele->measure() / (3.0 * ele->material->size);
	loc[0][0] =  val;
	loc[1][1] =  val;
	loc[0][1] = - val / 2.0;
	loc[1][0] = - val / 2.0;

	ele->bas_alfa[0] = ele->bas_alfa[1] = -1.0 / ele->measure();
	ele->bas_beta[0] = 1.0;
	ele->bas_beta[1] = 0.0;
}
//=============================================================================
// CALCULATE LOCAL MARIX FOR SIMPLEX ELEMENT
//=============================================================================
void local_matrix_tetrahedron( ElementIter ele )
{
	double alfa[ 4 ];
    double beta[ 4 ];       //  | Parametrs of basis functions
	double gama[ 4 ];       //  |
	double delta[ 4 ];      // /
	double poly[ 10 ];      // Koeficients of polynomial
	int i, j;               // Loops' counters
	SmallMtx4 loc=(SmallMtx4)(ele->loc);

	basis_functions_tetrahedron( ele, alfa, beta, gama, delta );
	for( i = 0; i < 4; i++ )
		for( j = i; j < 4; j++ ) {
			calc_polynom_tetrahedron(
					alfa[ i ], beta[ i ], gama[ i ],
					alfa[ j ], beta[ j ], gama[ j ],
					(SmallMtx3)(ele->a), poly );
			loc[i][j] = delta[i] * delta[j] * polynom_integral_tetrahedron( ele, poly );
			loc[j][i] = loc[i][j];
		}
        for(i = 0; i < 4; i++) {
                ele->bas_alfa[ i ] = alfa[ i ];
                ele->bas_beta[ i ] = beta[ i ];
                ele->bas_gama[ i ] = gama[ i ];
                ele->bas_delta[ i ] = delta[ i ];
        }
}
//=============================================================================
//
//=============================================================================
void basis_functions_tetrahedron( ElementIter ele, double alfa[], double beta[],double gama[], double delta[] )
{
	int li;
	struct Side *pSid;

	for( li = 0; li < 4; li++ ) {
		alfa[ li ] = ele->node[ li ]->getX();
		beta[ li ] = ele->node[ li ]->getY();
		gama[ li ] = ele->node[ li ]->getZ();
		pSid = ele->side[ li ];
		delta[ li ] = 1.0 / ( pSid->metric() *
			( pSid->normal[ 0 ] * pSid->centre()[ 0 ] +
			  pSid->normal[ 1 ] * pSid->centre()[ 1 ] +
			  pSid->normal[ 2 ] * pSid->centre()[ 2 ] -
			  pSid->normal[ 0 ] * alfa[ li ] -
			  pSid->normal[ 1 ] * beta[ li ] -
			  pSid->normal[ 2 ] * gama[ li ] ) );
	}
}
//=============================================================================
// CALCULATE POLYNOM OF SCALAR PRODUCT
//=============================================================================
void calc_polynom_tetrahedron( double al_i, double be_i, double ga_i,
                               double al_j, double be_j, double ga_j,
                               double a[ 3 ][ 3 ], double poly[] )
{
	// Constant term
        poly[ 0 ] =   a[ 0 ][ 0 ] * al_i * al_j +
		      a[ 0 ][ 1 ] * be_i * al_j +
		      a[ 0 ][ 2 ] * ga_i * al_j +
		      a[ 1 ][ 0 ] * al_i * be_j +
		      a[ 1 ][ 1 ] * be_i * be_j +
		      a[ 1 ][ 2 ] * ga_i * be_j +
		      a[ 2 ][ 0 ] * al_i * ga_j +
		      a[ 2 ][ 1 ] * be_i * ga_j +
		      a[ 2 ][ 2 ] * ga_i * ga_j;
	// Term with x
	poly[ 1 ] = ( a[ 0 ][ 0 ] * al_i +
		      a[ 0 ][ 1 ] * be_i +
		      a[ 0 ][ 2 ] * ga_i +
		      a[ 0 ][ 0 ] * al_j +
		      a[ 1 ][ 0 ] * be_j +
		      a[ 2 ][ 0 ] * ga_j ) * -1.0;
	// Term with y
	poly[ 2 ] = ( a[ 1 ][ 0 ] * al_i +
		      a[ 1 ][ 1 ] * be_i +
		      a[ 1 ][ 2 ] * ga_i +
		      a[ 0 ][ 1 ] * al_j +
		      a[ 1 ][ 1 ] * be_j +
		      a[ 2 ][ 1 ] * ga_j ) * -1.0;
	// Term with z
	poly[ 3 ] = ( a[ 2 ][ 0 ] * al_i +
		      a[ 2 ][ 1 ] * be_i +
		      a[ 2 ][ 2 ] * ga_i +
		      a[ 0 ][ 2 ] * al_j +
		      a[ 1 ][ 2 ] * be_j +
		      a[ 2 ][ 2 ] * ga_j ) * -1.0;
	// Term with xy
	poly[ 4 ] =   a[ 0 ][ 1 ] + a[ 1 ][ 0 ];
	// Term with xz
	poly[ 5 ] =   a[ 0 ][ 2 ] + a[ 2 ][ 0 ];
	// Term with yz
	poly[ 6 ] =   a[ 1 ][ 2 ] + a[ 2 ][ 1 ];
	// Term with x^2
	poly[ 7 ] =   a[ 0 ][ 0 ];
	// Term with y^2
	poly[ 8 ] =   a[ 1 ][ 1 ];
	// Term with z^2
	poly[ 9 ] =   a[ 2 ][ 2 ];
}
//=============================================================================
// CALCULATE INTEGRAL OF QUADRATICAL POLYNOM OVER TETRAHEDRON
//=============================================================================
double polynom_integral_tetrahedron( ElementIter ele, double poly[] )
{
	double rc;
	double v;
	double t[ 3 ];

	rc = 0.0;
	v = ele->measure();
	t[ 0 ] = ele->centre()[ 0 ];
	t[ 1 ] = ele->centre()[ 1 ];
	t[ 2 ] = ele->centre()[ 2 ];
	// Constant
	rc += poly[ 0 ] * v;
	// Term with x
	rc += poly[ 1 ] * v * t[ 0 ];
	// Term with y
	rc += poly[ 2 ] * v * t[ 1 ];
	// Term with z
	rc += poly[ 3 ] * v * t[ 2 ];
	// Term with xy
	rc += poly[ 4 ] * v * ( t[ 0 ] * t[ 1 ] + 1.0 / 20.0 * (
	      ( ele->node[ 0 ]->getX() - t[ 0 ] ) * ( ele->node[ 0 ]->getY() - t[ 1 ] ) +
	      ( ele->node[ 1 ]->getX() - t[ 0 ] ) * ( ele->node[ 1 ]->getY() - t[ 1 ] ) +
	      ( ele->node[ 2 ]->getX() - t[ 0 ] ) * ( ele->node[ 2 ]->getY() - t[ 1 ] ) +
	      ( ele->node[ 3 ]->getX() - t[ 0 ] ) * ( ele->node[ 3 ]->getY() - t[ 1 ] )
	      ) );
	// Term with xz
	rc += poly[ 5 ] * v * ( t[ 0 ] * t[ 2 ] + 1.0 / 20.0 * (
	      ( ele->node[ 0 ]->getX() - t[ 0 ] ) * ( ele->node[ 0 ]->getZ() - t[ 2 ] ) +
	      ( ele->node[ 1 ]->getX() - t[ 0 ] ) * ( ele->node[ 1 ]->getZ() - t[ 2 ] ) +
	      ( ele->node[ 2 ]->getX() - t[ 0 ] ) * ( ele->node[ 2 ]->getZ() - t[ 2 ] ) +
	      ( ele->node[ 3 ]->getX() - t[ 0 ] ) * ( ele->node[ 3 ]->getZ() - t[ 2 ] )
	      ) );
	// Term with yz
	rc += poly[ 6 ] * v * ( t[ 1 ] * t[ 2 ] + 1.0 / 20.0 * (
	      ( ele->node[ 0 ]->getY() - t[ 1 ] ) * ( ele->node[ 0 ]->getZ() - t[ 2 ] ) +
	      ( ele->node[ 1 ]->getY() - t[ 1 ] ) * ( ele->node[ 1 ]->getZ() - t[ 2 ] ) +
	      ( ele->node[ 2 ]->getY() - t[ 1 ] ) * ( ele->node[ 2 ]->getZ() - t[ 2 ] ) +
	      ( ele->node[ 3 ]->getY() - t[ 1 ] ) * ( ele->node[ 3 ]->getZ() - t[ 2 ] )
	      ) );
	// Term with x^2
	rc += poly[ 7 ] * v * ( t[ 0 ] * t[ 0 ] + 1.0 / 20.0 * (
	      ( ele->node[ 0 ]->getX() - t[ 0 ] ) * ( ele->node[ 0 ]->getX() - t[ 0 ] ) +
	      ( ele->node[ 1 ]->getX() - t[ 0 ] ) * ( ele->node[ 1 ]->getX() - t[ 0 ] ) +
	      ( ele->node[ 2 ]->getX() - t[ 0 ] ) * ( ele->node[ 2 ]->getX() - t[ 0 ] ) +
	      ( ele->node[ 3 ]->getX() - t[ 0 ] ) * ( ele->node[ 3 ]->getX() - t[ 0 ] )
	      ) );
	// Term with y^2
	rc += poly[ 8 ] * v * ( t[ 1 ] * t[ 1 ] + 1.0 / 20.0 * (
	      ( ele->node[ 0 ]->getY() - t[ 1 ] ) * ( ele->node[ 0 ]->getY() - t[ 1 ] ) +
	      ( ele->node[ 1 ]->getY() - t[ 1 ] ) * ( ele->node[ 1 ]->getY() - t[ 1 ] ) +
	      ( ele->node[ 2 ]->getY() - t[ 1 ] ) * ( ele->node[ 2 ]->getY() - t[ 1 ] ) +
	      ( ele->node[ 3 ]->getY() - t[ 1 ] ) * ( ele->node[ 3 ]->getY() - t[ 1 ] )
	      ) );
	// Term with z^2
	rc += poly[ 9 ] * v * ( t[ 2 ] * t[ 2 ] + 1.0 / 20.0 * (
	      ( ele->node[ 0 ]->getZ() - t[ 2 ] ) * ( ele->node[ 0 ]->getZ() - t[ 2 ] ) +
	      ( ele->node[ 1 ]->getZ() - t[ 2 ] ) * ( ele->node[ 1 ]->getZ() - t[ 2 ] ) +
	      ( ele->node[ 2 ]->getZ() - t[ 2 ] ) * ( ele->node[ 2 ]->getZ() - t[ 2 ] ) +
	      ( ele->node[ 3 ]->getZ() - t[ 2 ] ) * ( ele->node[ 3 ]->getZ() - t[ 2 ] )
	      ) );
	return rc;
}
//-----------------------------------------------------------------------------
// vim: set cindent:
