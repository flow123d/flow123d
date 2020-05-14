/*
 * fe_values_test.cpp
 *
 *  Created on: Sep 9, 2012
 *      Author: jb
 */

#define FEAL_OVERRIDE_ASSERTS

#include <flow_gtest.hh>
#include <cmath>
#include "arma_expect.hh"
#include "armadillo"
#include "system/armadillo_tools.hh"
#include "system/sys_profiler.hh"
#include "quadrature/quadrature_lib.hh"
#include "fem/fe_p.hh"
#include "fem/fe_values.hh"
#include "fem/mapping_p1.hh"
#include "mesh/mesh.h"
#include "mesh/elements.h"
#include "mesh/region.hh"
#include "mesh/accessors.hh"

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
double integrate(ElementAccessor<3> &ele) {
    FE_P_disc<dim> fe(0);
    QGauss quad( dim, 2 );
    MappingP1<dim,3> map;
    FEValues<3> fe_values(quad, fe, update_JxW_values | update_quadrature_points);
    
    fe_values.reinit( ele );
    
    double sum = 0.0;
    for(unsigned int i_point=0; i_point < fe_values.n_points(); i_point++) {
        sum += func<dim>( quad.point<dim>(i_point) ) * fe_values.JxW(i_point);
    }
    return sum;
}


TEST(FeValues, test_all) {
  // integrate a polynomial defined on the ref. element over an arbitrary element
    {
        // 1d case interval (1,3)   det(jac) = 2
    	Mesh mesh;
    	mesh.init_node_vector(2);
    	mesh.add_node(0, arma::vec3("1 0 0"));
    	mesh.add_node(1, arma::vec3("3 0 0"));
    	std::vector<unsigned int> node_ids = {0, 1};
    	mesh.init_element_vector(1);
    	mesh.add_element(0, 1, 1, 0, node_ids);
    	ElementAccessor<3> elm_acc = mesh.element_accessor(0);

        EXPECT_DOUBLE_EQ( 2.5 * 2, integrate<1>( elm_acc ) );
        
        // projection methods
        MappingP1<1,3> mapping;
        arma::mat::fixed<3, 2> map = mapping.element_map(elm_acc);
        EXPECT_ARMA_EQ( arma::mat("1 3; 0 0; 0 0"), map);
        EXPECT_ARMA_EQ( arma::vec("0.5 0.5"), mapping.project_real_to_unit( arma::vec("2.0 0.0 0.0"), map ) );
    }

    {
        // 2d case: triangle (0,1) (2,0) (3,4) surface = 3*4 - 1*2/2 - 1*4/4 - 3*3/2 = 9/2, det(jac) = 9
    	Mesh mesh;
    	mesh.init_node_vector(3);
    	mesh.add_node(0, arma::vec3("0 1 0"));
    	mesh.add_node(1, arma::vec3("2 0 0"));
    	mesh.add_node(2, arma::vec3("3 4 0"));
    	std::vector<unsigned int> node_ids = {0, 1, 2};
    	mesh.init_element_vector(1);
    	mesh.add_element(1, 2, 1, 0, node_ids);
    	ElementAccessor<3> elm_acc = mesh.element_accessor(0);

        EXPECT_DOUBLE_EQ( 19.0 / 12.0 * 9.0 , integrate<2>( elm_acc ) );

        // projection methods
        MappingP1<2,3> mapping;
        arma::mat::fixed<3, 3> map = mapping.element_map(elm_acc);
        EXPECT_ARMA_EQ( arma::mat("0 2 3; 1 0 4; 0 0 0"), map);
        EXPECT_ARMA_EQ( arma::vec("0.6 0.2 0.2"), mapping.project_real_to_unit( arma::vec("1.0 1.4 0.0"), map ) );
    }

}


class TestElementMapping {
public:
    TestElementMapping(std::vector<string> nodes_str)
    {
        unsigned int i=0;
    	mesh_.init_node_vector(4);
    	mesh_.init_element_vector(1);
    	for(auto str : nodes_str) mesh_.add_node(i++, arma::vec3(str));
    	std::vector<unsigned int> node_ids = {0, 1, 2, 3};
    	mesh_.add_element(1, 3, 1, 0, node_ids);
    }

    ElementAccessor<3> elem_accessor()
	{
    	return mesh_.element_accessor(0);
	}

    Mesh mesh_;
};


TEST(ElementMapping, element_map) {
    Profiler::instance();
    armadillo_setup();
    MappingP1<3,3> mapping;

    {
        TestElementMapping ele_mapping({ "0 0 0", "1 0 0", "0 1 0", "0 0 1"});
        arma::mat::fixed<3, 4> map = mapping.element_map(ele_mapping.elem_accessor());
        EXPECT_ARMA_EQ( arma::mat("0 1 0 0; 0 0 1 0; 0 0 0 1"), map);
        EXPECT_ARMA_EQ( arma::vec("0.4 0.1 0.2 0.3"), mapping.project_real_to_unit( arma::vec3("0.1 0.2 0.3"), map ) );
        EXPECT_ARMA_EQ( arma::vec("-0.5 0.5 0.5 0.5"), mapping.project_real_to_unit( arma::vec3("0.5 0.5 0.5"), map ) );
    }

    {
        // trnaslated
        TestElementMapping ele_mapping({ "1 2 3", "2 2 3", "1 3 3", "1 2 4"});
        arma::mat::fixed<3, 4> map = mapping.element_map(ele_mapping.elem_accessor());
        EXPECT_ARMA_EQ( arma::mat("1 2 1 1; 2 2 3 2; 3 3 3 4"), map);
        EXPECT_ARMA_EQ( arma::vec("0.4 0.1 0.2 0.3"), mapping.project_real_to_unit( arma::vec3("1.1 2.2 3.3"), map ) );
        EXPECT_ARMA_EQ( arma::vec("-0.5 0.5 0.5 0.5"), mapping.project_real_to_unit( arma::vec3("1.5 2.5 3.5"), map ) );
    }

    {
        // simplest cube element 7
        TestElementMapping ele_mapping({ "-1 -1 1", "1 1 -1", "-1 -1 -1", "1 -1 -1"});
        arma::mat::fixed<3, 4> map = mapping.element_map(ele_mapping.elem_accessor());
        EXPECT_ARMA_EQ( arma::mat("-1 1 -1 1; -1 1 -1 -1; 1 -1 -1 -1"), map);
        EXPECT_ARMA_EQ( arma::vec("0.25 0.25 0.25 0.25"), mapping.project_real_to_unit( arma::vec3("0 -0.5 -0.5"), map ) );
        //EXPECT_ARMA_EQ( arma::vec("0.1 0.2 0.3 0.4"), mapping.project_real_to_unit( arma::vec3("0.1 0.2 0.3"), map ) );
    }
}
