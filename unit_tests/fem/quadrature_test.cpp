/*
 * python_function_test.cpp
 *
 *  Created on: Aug 30, 2012
 *      Author: jb
 */



#include <flow_gtest.hh>
#include <cmath>
#include "quadrature/quadrature_lib.hh"
#include "quadrature/qmidpoint.hh"
#include "arma_expect.hh"

#define INTEGRATE( _func_ ) for( unsigned int i=0; i < quad.size(); i++) sum +=  _func_( quad.point(i) ) * quad.weight(i);

double test_1_1d( const arma::vec::fixed<1> & p) {
    return 3 * p[0] + 1.0;
}

double test_2_1d( const arma::vec::fixed<1> & p) {
    return 3 * p[0] * p[0] + p[0] + 1.0;
}


TEST(Quadrature, test_1d) {

    {
    QGauss<1> quad( 1 ); // should integrate P1 exactly
    EXPECT_EQ(1, quad.size());
    double sum =0.0;
    INTEGRATE(test_1_1d);
    EXPECT_DOUBLE_EQ(5.0/2.0, sum); // 3 * 1/2 + 1
    }

    {
    QGauss<1> quad( 2 ); // should integrate P2 exactly
    EXPECT_EQ(2, quad.size());
    double sum =0.0;
    INTEGRATE(test_2_1d);
    EXPECT_DOUBLE_EQ(2.5, sum); // 3 * 1/3 + 1/2 + 1
    }
}





double test_1_2d( const arma::vec::fixed<2> & p) {
    return 3 * p[1] + 2 * p[0] + 1.0;
}

double test_2_2d( const arma::vec::fixed<2> & p) {
    return 3 * p[0] * p[0] + p[0] + 6 * p[1] * p[1] + p[1] + 1.0;
}

TEST(Quadrature, test_2d) {

    {
    QGauss<2> quad( 1 ); // should integrate P1 exactly
    EXPECT_EQ(1, quad.size());
    double sum =0.0;
    INTEGRATE(test_1_2d);
    EXPECT_DOUBLE_EQ(8.0/6.0, sum); // 3 * 1/6 + 2 * 1/6 + 1/2 = 8/6
    }

    {
    QGauss<2> quad( 2 ); // should integrate P2 exactly
    EXPECT_EQ(3, quad.size());
    double sum =0.0;
    INTEGRATE(test_2_2d);
    EXPECT_DOUBLE_EQ(19.0 / 12.0 , sum); // 3 * 1/12 + 1/6 + 6 * 1/12 + 1/6 + 1/2 = 19/12
    }
}


TEST(Quadrature, midpoint){
    QMidpoint quad(25);
    EXPECT_EQ(25, quad.size());
    
    double sum =0.0;
    INTEGRATE(test_1_1d);   // should integrate P1 exactly
    EXPECT_DOUBLE_EQ(5.0/2.0, sum); // 3 * 1/2 + 1
}


TEST(Quadrature, side_projections){
    {
    QMidpoint qm(4);
    
    Quadrature<2> quad(qm,0,0);
    EXPECT_ARMA_EQ(arma::vec("0.125 0"), quad.point(0));
    EXPECT_ARMA_EQ(arma::vec("0.375 0"), quad.point(1));
    EXPECT_ARMA_EQ(arma::vec("0.625 0"), quad.point(2));
    EXPECT_ARMA_EQ(arma::vec("0.875 0"), quad.point(3));
    quad = Quadrature<2>(qm,1,0);
    EXPECT_ARMA_EQ(arma::vec("0 0.125"), quad.point(0));
    EXPECT_ARMA_EQ(arma::vec("0 0.625"), quad.point(2));
    quad = Quadrature<2>(qm,2,0);
    EXPECT_ARMA_EQ(arma::vec("0.875 0.125"), quad.point(0));
    EXPECT_ARMA_EQ(arma::vec("0.375 0.625"), quad.point(2));
    quad = Quadrature<2>(qm,2,1);
    EXPECT_ARMA_EQ(arma::vec("0.125 0.875"), quad.point(0));
    EXPECT_ARMA_EQ(arma::vec("0.625 0.375"), quad.point(2));
    }
    
    {
    QGauss<2> qg(1);
    Quadrature<3> quad(qg, 0, 0);
    EXPECT_ARMA_EQ(arma::vec({1./3, 1./3, 0}), quad.point(0));
    quad = Quadrature<3>(qg, 1, 0);
    EXPECT_ARMA_EQ(arma::vec({1./3, 0, 1./3}), quad.point(0));
    quad = Quadrature<3>(qg, 2, 0);
    EXPECT_ARMA_EQ(arma::vec({0, 1./3, 1./3}), quad.point(0));
    quad = Quadrature<3>(qg, 3, 0);
    EXPECT_ARMA_EQ(arma::vec({1./3, 1./3, 1./3}), quad.point(0));
    }
}

