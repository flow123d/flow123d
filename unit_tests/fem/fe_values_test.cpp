/*
 * fe_values_test.cpp
 *
 *  Created on: Sep 9, 2012
 *      Author: jb
 */

#include <flow_gtest.hh>
#include <cmath>
#include "quadrature/quadrature_lib.hh"
#include "fem/fe_p.hh"
#include "fem/fe_values.hh"
#include "fem/mapping_p1.hh"
#include "mesh/mesh.h"

#define INTEGRATE( _func_ ) for( unsigned int i=0; i < quad.size(); i++) sum +=  _func_( quad.point(i) ) * quad.weight(i);

double test_1_1d( const arma::vec::fixed<1> & p) {
    return 3 * p[0] + 1.0;
}

template <int dim>
double func( const arma::vec::fixed<dim> & p) {
    if (dim == 1) {
        return 3 * p[0] * p[0] + p[0] + 1.0;
    } else {
        return 3 * p[0] * p[0] + p[0] + 6 * p[1] * p[1] + p[1] + 1.0;
    }
}


template <int dim>
double integrate(ElementFullIter &ele) {
    FE_P_disc<0,dim,3> fe;
    QGauss<dim> quad( 2 );
    MappingP1<dim,3> map;
    FEValues<dim,3> fe_values(map, quad,   fe, update_JxW_values | update_quadrature_points);

    fe_values.reinit( ele );

    double sum = 0.0;
    for(unsigned int i_point=0; i_point < fe_values.n_points(); i_point++) {
        sum += func<dim>( quad.point(i_point) ) * fe_values.JxW(i_point);
    }
    return sum;
}


TEST(FeValues, test_all) {
  // integrate a polynomial defined on the ref. element over an arbitrary element
    {
        // 1d case interval (1,3)   det(jac) = 2
        NodeVector nodes(2);
        nodes.add_item(0);
        nodes[0].point()[0] = 1.0;
        nodes[0].point()[1] = 0.0;
        nodes[0].point()[2] = 0.0;

        nodes.add_item(1);
        nodes[1].point()[0] = 3.0;
        nodes[1].point()[1] = 0.0;
        nodes[1].point()[2] = 0.0;

        ElementVector el_vec(1);
        el_vec.add_item(0);

        RegionIdx reg;
        Element ele(1, NULL, reg);      //NULL - mesh pointer, empty RegionIdx

        ele.node = new Node * [ele.n_nodes()];
        for(int i =0; i < 2; i++) ele.node[i] = nodes(i);
        el_vec[0] = ele; // dangerous since Element has no deep copy constructor.


        ElementFullIter it( el_vec(0) );
        EXPECT_DOUBLE_EQ( 2.5 * 2, integrate<1>( it ) );

    }

    {
        // 2d case: triangle (0,1) (2,0) (3,4) surface = 3*4 - 1*2/2 - 1*4/4 - 3*3/2 = 9/2, det(jac) = 9
        NodeVector nodes(3);
        nodes.add_item(0);
        nodes[0].point()[0] = 0.0;
        nodes[0].point()[1] = 1.0;
        nodes[0].point()[2] = 0.0;

        nodes.add_item(1);
        nodes[1].point()[0] = 2.0;
        nodes[1].point()[1] = 0.0;
        nodes[1].point()[2] = 0.0;

        nodes.add_item(2);
        nodes[2].point()[0] = 3.0;
        nodes[2].point()[1] = 4.0;
        nodes[2].point()[2] = 0.0;

        ElementVector el_vec(1);
        el_vec.add_item(0);

        RegionIdx reg;
        Element ele(2, NULL, reg);      //NULL - mesh pointer, empty RegionIdx

        ele.node = new Node * [ele.n_nodes()];
        for(int i =0; i < 3; i++) ele.node[i] = nodes(i);
        el_vec[0] = ele; // dangerous since Element has no deep copy constructor.

        ElementFullIter it( el_vec(0) );
        EXPECT_DOUBLE_EQ( 19.0 / 12.0 * 9.0 , integrate<2>( it ) );
    }

}


