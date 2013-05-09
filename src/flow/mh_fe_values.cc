/*
 * mh_fe_values.cc
 *
 *  Created on: Jul 2, 2012
 *      Author: jb
 */


#include "system/system.hh"
#include "system/math_fce.h"
#include "mesh/mesh.h"

#include "mh_fe_values.hh"


MHFEValues::MHFEValues() {
    // allocate temp. arrays
    loc_matrix_ = new double [4*4 + 4*4 + 4*4];
    inv_loc_matrix_ = loc_matrix_ + 4*4;
    bas_alfa = inv_loc_matrix_ + 4*4;
    bas_beta = bas_alfa + 4;
    bas_gama = bas_beta +4;
    bas_delta = bas_gama +4;
}



MHFEValues::~MHFEValues() {
    delete [] loc_matrix_;
}



void MHFEValues::update(ElementFullIter ele, FieldType &anisotropy, FieldType_Scalar &cross_section, FieldType_Scalar &conductivity) {

    ASSERT(!( ele == NULL ),"NULL as argument of function local_matrix()\n");

    double scale = 1/ conductivity.value( ele->centre(), ele->element_accessor() ) / cross_section.value( ele->centre(), ele->element_accessor() );
    //DBGMSG("scale: %g\n", scale);
    switch( ele->dim() ) {
        case 1:
            local_matrix_line( ele, anisotropy , scale);
            break;
        case 2:
            local_matrix_triangle( ele, anisotropy, scale);
            break;
        case 3:
            local_matrix_tetrahedron( ele, anisotropy, scale );
            break;
    }

    // matrix inverse

    double det = MatrixInverse(loc_matrix_, inv_loc_matrix_, ele->n_sides() );
    if (fabs(det) < NUM_ZERO) {
        xprintf(Warn,"Singular local matrix of the element %d\n",ele.id());
        PrintSmallMatrix(loc_matrix_, ele->n_sides());
        xprintf(Err,"det: %30.18e \n",det);
    }
}



double * MHFEValues::local_matrix() {
    return loc_matrix_;
}



double * MHFEValues::inv_local_matrix() {
    return inv_loc_matrix_;
}



arma::vec3 MHFEValues::RT0_value(ElementFullIter ele, arma::vec3 point, unsigned int face) {
    switch( ele->dim() ) {
    case 1:
        {

            arma::vec3 line_vec = ele->node[1]->point() - ele->node[0]->point();
            if (face == 0) {
                return - arma::norm( point - ele->node[1]->point(), 2) * line_vec / arma::dot( line_vec, line_vec) ;
            } else {
                return arma::norm( point - ele->node[0]->point(), 2) * line_vec / arma::dot( line_vec, line_vec) ;
            }
        }
    case 2:
        {
            // make rotated coordinate system with triangle in plane XY, origin in A and axes X == AB
            arma::vec3 ex(ele->node[1]->point() - ele->node[0]->point());
            ex /= norm(ex,2);

            arma::vec3 ac(ele->node[2]->point() - ele->node[0]->point());
            arma::vec3 ez = cross(ex, ac);
            ez /= norm(ez,2);

            arma::vec3 ey = cross(ez,ex);
            ey /= norm(ey, 2);

            // compute point in new coordinate system
            arma::vec3 u = point - ele->node[0]->point();

            // compute vector value from the base function
            return bas_gama[ face ] * (
                         (dot(u, ex) - bas_alfa[ face ] ) * ex
                        +(dot(u, ey) - bas_beta[ face ] ) * ey
                   );
        }
    case 3:
        {
            arma::vec3 RT0_Y;
            RT0_Y[0] = bas_alfa[ face ];
            RT0_Y[1] = bas_beta[ face ];
            RT0_Y[2] = bas_gama[ face ];

            return bas_delta[ face ]*(point - RT0_Y);
        }
    }

    return arma::vec3();
}





//=============================================================================
// CALCULATE LOCAL MARIX FOR LINEAR ELEMENT
//=============================================================================
void MHFEValues::local_matrix_line(ElementFullIter ele, FieldType &anisotropy, double scale )
{
    double    val;
    SmallMtx2 loc=(SmallMtx2)(loc_matrix_);
    
    //getting vector on the line and normalizing it
    //it is the transformation matrix from 3D to 1D line of the 1D element
    arma::vec3 line_vec = ele->node[1]->point() - ele->node[0]->point();
    line_vec /= arma::norm(line_vec,2);
    
    //transforming the conductivity in 3D to resistivity in 1D
    //computing v_transpose * K_inverse * v
    val = scale * arma::dot(line_vec,
                        (anisotropy.value(ele->centre(), ele->element_accessor() )).i() * line_vec
                    ) * ele->measure() / 3.0;
              
    loc[0][0] =  val;
    loc[1][1] =  val;
    loc[0][1] = - val / 2.0;
    loc[1][0] = - val / 2.0;

    bas_alfa[0] = bas_alfa[1] = -1.0 / ele->measure();
    bas_beta[0] = 1.0;
    bas_beta[1] = 0.0;
}


//=============================================================================
// CALCULATE LOCAL MATRIX FOR TRIANGULAR ELEMENT
//=============================================================================
void MHFEValues::local_matrix_triangle( ElementFullIter ele, FieldType &anisotropy, double scale )
{
    double midpoint[ 3 ][ 2 ]; // Midpoints of element's sides
    double alfa[ 3 ];   //
    double beta[ 3 ];       //  |- Parametrs of basis functions
    double gama[ 3 ];       // /
    double poly[ 6 ];       // Koeficients of polynomial
    double p[ 3 ];          // Polynomial's values in midpoints
    int i, j;               // Loops' counters
    double nod_coor[ 3 ][ 2 ]; // Coordinates of nodes;
    SmallMtx3 loc=(SmallMtx3)(loc_matrix_);
    
    // make rotated coordinate system with triangle in plane XY, origin in A and axes X == AB
    arma::vec3 ex(ele->node[1]->point() - ele->node[0]->point());
    ex /= arma::norm(ex,2);

    arma::vec3 ac(ele->node[2]->point() - ele->node[0]->point());
    arma::vec3 ez = arma::cross(ex, ac);
    ez /= norm(ez,2);

    arma::vec3 ey = arma::cross(ez,ex);
    ey /= arma::norm(ey, 2);
            
    //transformation matrix from 3D to 2D plane of the 2D element
    arma::mat r(3,2);
    r.col(0) = ex;
    r.col(1) = ey;
    
    //transforming 3D conductivity tensor to 2D resistance tensor
    arma::mat resistance_tensor = r.t() * ((anisotropy.value(ele->centre(), ele->element_accessor() )).i() * r);

    node_coordinates_triangle( ele, nod_coor );
    side_midpoint_triangle( nod_coor, midpoint );
    basis_functions_triangle( nod_coor, alfa, beta, gama );
    for( i = 0; i < 3; i++ )
        for( j = i; j < 3; j++ ) {
            calc_polynom_triangle( alfa[ i ], beta[ i ],
                               alfa[ j ], beta[ j ],
                             resistance_tensor, poly );
            p[ 0 ] = polynom_value_triangle( poly, midpoint[ 0 ] );
            p[ 1 ] = polynom_value_triangle( poly, midpoint[ 1 ] );
            p[ 2 ] = polynom_value_triangle( poly, midpoint[ 2 ] );
            loc[i][j] = scale * gama[i] * gama[j] * ele->measure() *
                    ( p[0] + p[1] + p[2] ) / 3.0;
            loc[j][i] = loc[i][j];
        }
        for(i = 0; i < 3; i++) {
                bas_alfa[ i ] = alfa[ i ];
                bas_beta[ i ] = beta[ i ];
                bas_gama[ i ] = gama[ i ];
        }
    //check_local( ele );
}


/*
 * Computes coordinates of vertices of the triangle  in the local orthogonal system
 * that has normalized vector  V0 to V1 as the first vector of the basis.
 */
void MHFEValues::node_coordinates_triangle( ElementFullIter ele, double nod[ 3 ][ 2 ] )
{
    double u[ 3 ], v[ 3 ];  // vectors of sides 0 and 2
    double t[ 3 ];      // u normalized (tangenta)
    double n[ 3 ];      // normal
    double b[ 3 ];      // binormal

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


/**
 * Computes coordinates of the side midpoints in the local orthogonal coordinate system.
 */
void MHFEValues::side_midpoint_triangle( double nod[ 3 ][ 2 ], double midpoint[ 3 ][ 2 ] )
{
    midpoint[ 0 ][ 0 ] = ( nod[ 0 ][ 0 ] + nod[ 1 ][ 0 ] ) / 2.0;
    midpoint[ 0 ][ 1 ] = ( nod[ 0 ][ 1 ] + nod[ 1 ][ 1 ] ) / 2.0;
    midpoint[ 1 ][ 0 ] = ( nod[ 1 ][ 0 ] + nod[ 2 ][ 0 ] ) / 2.0;
    midpoint[ 1 ][ 1 ] = ( nod[ 1 ][ 1 ] + nod[ 2 ][ 1 ] ) / 2.0;
    midpoint[ 2 ][ 0 ] = ( nod[ 2 ][ 0 ] + nod[ 0 ][ 0 ] ) / 2.0;
    midpoint[ 2 ][ 1 ] = ( nod[ 2 ][ 1 ] + nod[ 0 ][ 1 ] ) / 2.0;
}



/*
 * Computes coefficients of the RT basis functions.
 */
void MHFEValues::basis_functions_triangle( double nod[ 3 ][ 2 ], double alfa[],
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
void MHFEValues::bas_func_0_triangle( double x0, double y0,
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
void MHFEValues::calc_polynom_triangle( double al_i, double be_i, double al_j, double be_j,
                        arma::mat::fixed<2,2> a, double poly[] )
{
        poly[ 0 ] =   a( 0,0 ) * al_i * al_j +
                      a( 0,1 ) * be_i * al_j +
                      a( 1,0 ) * al_i * be_j +
                      a( 1,1 ) * be_i * be_j;
                      
        poly[ 1 ] = ( a( 0,0 ) * al_i +
                      a( 0,0 ) * al_j +
                      a( 0,1 ) * be_i +
                      a( 1,0 ) * be_j ) * -1.0;
                      
        poly[ 2 ] = ( a( 1,1 ) * be_i +
                      a( 1,1 ) * be_j +
                      a( 1,0 ) * al_i +
                      a( 0,1) * al_j ) * -1.0;
        poly[ 3 ] =   a( 0,0 );
        poly[ 4 ] =   a( 0,1 ) + a( 1,0 );
        poly[ 5 ] =   a( 1,1 );
}
//=============================================================================
// CALCULATE VALUE OF POLYNOM IN GIVEN POINT
//=============================================================================
double MHFEValues::polynom_value_triangle( double poly[], double point[] )
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
// CALCULATE LOCAL MARIX FOR SIMPLEX ELEMENT
//=============================================================================
void MHFEValues::local_matrix_tetrahedron( ElementFullIter ele, FieldType &anisotropy, double scale )
{
    double alfa[ 4 ];
    double beta[ 4 ];       //  | Parametrs of basis functions
    double gama[ 4 ];       //  |
    double delta[ 4 ];      // /
    double poly[ 10 ];      // Koeficients of polynomial
    int i, j;               // Loops' counters
    SmallMtx4 loc=(SmallMtx4)(loc_matrix_);
    
    //transforming 3D conductivity tensor to 3D resistance tensor
    arma::mat resistance_tensor = (anisotropy.value(ele->centre(), ele->element_accessor() )).i();
        
    //OBSOLETE
    //SmallMtx3 resistance_tensor = (SmallMtx3)(ele->material->hydrodynamic_resistence);
    //compares old and new resistance tensor 
    //DBGMSG("my: %f %f %f %f\t orig: %f %f %f %f\n", 
    //        my_resistance_tensor(0,0), my_resistance_tensor(0,1), my_resistance_tensor(1,0), my_resistance_tensor(1,1), 
    //       resistance_tensor[0][0], resistance_tensor[0][1], resistance_tensor[1][0], resistance_tensor[1][1] );

    basis_functions_tetrahedron( ele, alfa, beta, gama, delta);
    for( i = 0; i < 4; i++ )
        for( j = i; j < 4; j++ ) {
            calc_polynom_tetrahedron(
                    alfa[ i ], beta[ i ], gama[ i ],
                    alfa[ j ], beta[ j ], gama[ j ],
                    resistance_tensor, poly );
            loc[i][j] = scale * delta[i] * delta[j] * polynom_integral_tetrahedron( ele, poly );
            loc[j][i] = loc[i][j];
        }
        for(i = 0; i < 4; i++) {
                bas_alfa[ i ] = alfa[ i ];
                bas_beta[ i ] = beta[ i ];
                bas_gama[ i ] = gama[ i ];
                bas_delta[ i ] = delta[ i ];
        }
}
//=============================================================================
//
//=============================================================================
void MHFEValues::basis_functions_tetrahedron( ElementFullIter ele, double alfa[], double beta[],double gama[], double delta[])
{
    int li;
    SideIter pSid;

    for( li = 0; li < 4; li++ ) {
        alfa[ li ] = ele->node[ li ]->getX();
        beta[ li ] = ele->node[ li ]->getY();
        gama[ li ] = ele->node[ li ]->getZ();
        pSid = ele->side( li );
        delta[ li ] = 1.0 / ( pSid->measure() * 
            ( pSid->normal()[ 0 ] * pSid->centre()[ 0 ] +
              pSid->normal()[ 1 ] * pSid->centre()[ 1 ] +
              pSid->normal()[ 2 ] * pSid->centre()[ 2 ] -
              pSid->normal()[ 0 ] * alfa[ li ] -
              pSid->normal()[ 1 ] * beta[ li ] -
              pSid->normal()[ 2 ] * gama[ li ] ) );
    }
}
//=============================================================================
// CALCULATE POLYNOM OF SCALAR PRODUCT
//=============================================================================
void MHFEValues::calc_polynom_tetrahedron( double al_i, double be_i, double ga_i,
                               double al_j, double be_j, double ga_j,
                               arma::mat::fixed<3,3> a, double poly[] )
{
    // Constant term
  poly[ 0 ] =   a( 0,0 ) * al_i * al_j +
                a( 0,1 ) * be_i * al_j +
                a( 0,2 ) * ga_i * al_j +
                a( 1,0 ) * al_i * be_j +
                a( 1,1 ) * be_i * be_j +
                a( 1,2 ) * ga_i * be_j +
                a( 2,0 ) * al_i * ga_j +
                a( 2,1 ) * be_i * ga_j +
                a( 2,2 ) * ga_i * ga_j;
    // Term with x
  poly[ 1 ] = ( a( 0,0 ) * al_i +
                a( 0,1 ) * be_i +
                a( 0,2 ) * ga_i +
                a( 0,0 ) * al_j +
                a( 1,0 ) * be_j +
                a( 2,0 ) * ga_j ) * -1.0;
    // Term with y
  poly[ 2 ] = ( a( 1,0 ) * al_i +
                a( 1,1 ) * be_i +
                a( 1,2 ) * ga_i +
                a( 0,1 ) * al_j +
                a( 1,1 ) * be_j +
                a( 2,1 ) * ga_j ) * -1.0;
    // Term with z
  poly[ 3 ] = ( a( 2,0 ) * al_i +
                a( 2,1 ) * be_i +
                a( 2,2 ) * ga_i +
                a( 0,2 ) * al_j +
                a( 1,2 ) * be_j +
                a( 2,2 ) * ga_j ) * -1.0;
    // Term with xy
    poly[ 4 ] =   a( 0,1 ) + a( 1,0 );
    // Term with xz
    poly[ 5 ] =   a( 0,2 ) + a( 2,0 );
    // Term with yz
    poly[ 6 ] =   a( 1,2 ) + a( 2,1 );
    // Term with x^2
    poly[ 7 ] =   a( 0,0 );
    // Term with y^2
    poly[ 8 ] =   a( 1,1 );
    // Term with z^2
    poly[ 9 ] =   a( 2,2 );
}
//=============================================================================
// CALCULATE INTEGRAL OF QUADRATICAL POLYNOM OVER TETRAHEDRON
//=============================================================================
double MHFEValues::polynom_integral_tetrahedron( ElementFullIter ele, double poly[] )
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
