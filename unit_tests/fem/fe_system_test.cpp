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
#include "mesh/element_impls.hh"
#include "mesh/region.hh"
#include "fem/fe_system.hh"

#define INTEGRATE( _func_ ) for( unsigned int i=0; i < quad.size(); i++) sum +=  _func_( quad.point(i) ) * quad.weight(i);


NodeVector make_nodes(const std::vector<string> &nodes_str)
{
  std::vector<arma::vec3> nodes;
  for(auto str : nodes_str) nodes.push_back( arma::vec3(str));
  
  NodeVector node_vector(nodes.size());
  unsigned int i=0;
  for (auto node : nodes)
  {
    node_vector.add_item(i);
    node_vector[i++] = Node(node[0], node[1], node[2]);
  }
  
  return node_vector;
}


ElementVector make_elements(NodeVector &node_vector, const std::vector<std::vector<unsigned int> > &node_idx)
{
  ElementVector el_vec(node_idx.size());
  
  unsigned int iel = 0;
  for (auto nodes : node_idx)
  {
    Element e;
    el_vec.add_item(iel);
    el_vec[iel].init(nodes.size()-1, nullptr, RegionIdx());
    unsigned int i=0;
    for(auto node : nodes)
      el_vec[iel].node[i++] = &node_vector[node];
    
    iel++;
  }
  
  return el_vec;
};





TEST(FeSystem, test_all) {
    // Test vector-valued FESystem using P1 element on tetrahedron.
    {
        NodeVector nodes = make_nodes({"0 0 0", "1 0 0", "0 1 0", "0 0 1"});
        ElementVector el_vec = make_elements(nodes, { { 0, 1, 2, 3 } });
        ElementFullIter ele( el_vec(0) );
        
        FESystem<3,3> fe_sys(new FE_P<1,3,3>, 3);
        MappingP1<3,3> map;
        Quadrature<3> q(nodes.size());
        for (unsigned int i=0; i<nodes.size(); i++)
          q.set_point(i, nodes[i].point());
        FEValues<3,3> fe_values(map, q, fe_sys, update_values | update_gradients);
        FEValuesExtractors::Vector vec(0);
        
        fe_values.reinit(ele);
        
        EXPECT_ARMA_EQ( arma::vec("1 0 0"), fe_values[vec].value(0,0));
        EXPECT_ARMA_EQ( arma::vec("0 1 0"), fe_values[vec].value(1,0));
        EXPECT_ARMA_EQ( arma::vec("0 0 1"), fe_values[vec].value(2,0));
        
        EXPECT_ARMA_EQ( arma::vec("1 0 0"), fe_values[vec].value(3,1));
        EXPECT_ARMA_EQ( arma::vec("0 1 0"), fe_values[vec].value(4,1));
        EXPECT_ARMA_EQ( arma::vec("0 0 1"), fe_values[vec].value(5,1));
        
        EXPECT_ARMA_EQ( arma::vec("1 0 0"), fe_values[vec].value(6,2));
        EXPECT_ARMA_EQ( arma::vec("0 1 0"), fe_values[vec].value(7,2));
        EXPECT_ARMA_EQ( arma::vec("0 0 1"), fe_values[vec].value(8,2));
        
        EXPECT_ARMA_EQ( arma::vec("1 0 0"), fe_values[vec].value(9,3));
        EXPECT_ARMA_EQ( arma::vec("0 1 0"), fe_values[vec].value(10,3));
        EXPECT_ARMA_EQ( arma::vec("0 0 1"), fe_values[vec].value(11,3));
        

    }

}



