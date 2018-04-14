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
#include "mesh/side_impl.hh"
#include "fem/fe_p.hh"
#include "fem/fe_rt.hh"
#include "fem/fe_values.hh"
#include "fem/mapping_p1.hh"
#include "mesh/mesh.h"
#include "mesh/elements.h"
#include "mesh/region.hh"
#include "mesh/accessors.hh"
#include "fem/fe_system.hh"



/*NodeVector make_nodes(const std::vector<string> &nodes_str)
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


vector<Element> make_elements(NodeVector &node_vector, const std::vector<std::vector<unsigned int> > &node_idx)
{
  vector<Element> el_vec(node_idx.size());
  
  unsigned int iel = 0;
  for (auto nodes : node_idx)
  {
    el_vec[iel].init(nodes.size()-1, iel, nullptr, RegionIdx());
    unsigned int i=0;
    for(auto node : nodes)
      el_vec[iel].node[i++] = &node_vector[node];
    
    iel++;
  }
  
  return el_vec;
}*/




class FESystemTest : public testing::Test {
public:
  FESystemTest()
    : q(4)
  {
  	mesh.add_node(0, arma::vec3("1 0 0"));
  	mesh.add_node(1, arma::vec3("0 1 0"));
  	mesh.add_node(2, arma::vec3("0 0 1"));
  	mesh.add_node(3, arma::vec3("0 0 0"));
  	std::vector<unsigned int> node_ids = {3, 0, 1, 2};
  	mesh.add_element(0, 3, 1, 0, node_ids);
  	ele = mesh.element_accessor(0);

  	for (unsigned int i=0; i<mesh.n_nodes(); i++)
      q.set_point(i, mesh.node_vector[i].point());
  }
  
protected:
  Mesh mesh;
  ElementAccessor<3> ele;
  MappingP1<3,3> map;
  Quadrature<3> q;

};









TEST_F(FESystemTest, test_vector) {
  // Test vector-valued FESystem using P1 element on tetrahedron.
  FESystem<3> fe_sys(std::make_shared<FE_P<3> >(1), FEVectorContravariant);
  FEValues<3,3> fe_values(map, q, fe_sys, update_values | update_gradients);
  
  fe_values.reinit(ele);
  
  auto vec_view = fe_values.vector_view(0);
  
  for (unsigned int k=0; k<q.size(); k++)
    for (unsigned int i=0; i<fe_sys.n_dofs(); i++)
      for (unsigned int c=0; c<3; c++)
      {
        // check values
        EXPECT_EQ( ((i%4==(k+1)%4) && (i/4==c))?1:0, vec_view.value(i,k)[c] );
        //check gradients
        arma::rowvec gr = vec_view.grad(i,k).row(c);
        if (i / 4 == c)
        { // gradient of nonzero component
          switch (i%4)
          {
            case 0:
              EXPECT_ARMA_EQ( arma::rowvec("-1 -1 -1"), gr );
              break;
            case 1:
              EXPECT_ARMA_EQ( arma::rowvec("1 0 0"), gr );
              break;
            case 2:
              EXPECT_ARMA_EQ( arma::rowvec("0 1 0"), gr );
              break;
            case 3:
              EXPECT_ARMA_EQ( arma::rowvec("0 0 1"), gr );
              break;
          }
        }
        else
          EXPECT_ARMA_EQ( arma::rowvec("0 0 0"), gr );
      }
}


TEST_F(FESystemTest, test_mixed_system) {
  // Test mixed-system FE using P0, P1^3 and RT0 elements on tetrahedron.
  // The basis functions are ordered first nodal and then element-supported,
  // hence the scalar constant function from P0 comes after the linear
  // functions from P1^3 and the RT0 functions are at the end.
  FESystem<3> fe_vec(std::make_shared<FE_P<3> >(1), FEVectorContravariant);
  FESystem<3> fe_sys({ std::make_shared<FE_P<3> >(0), std::make_shared<FESystem<3> >(fe_vec), std::make_shared<FE_RT0<3> >() });
  FEValues<3,3> fe_values(map, q, fe_sys, update_values | update_gradients);
  
  fe_values.reinit(ele);
  
  auto vec_view = fe_values.vector_view(0);
  auto scalar_view = fe_values.scalar_view(0);
  auto rt_view = fe_values.vector_view(1);
  
  for (unsigned int k=0; k<q.size(); k++)
  {
    // check values and gradients of P1^3 function
    for (unsigned int i=1; i<1+fe_vec.n_dofs(); i++)
      for (unsigned int c=0; c<3; c++)
      {
        // check values
        EXPECT_EQ( (((i-k)%4==2) && ((i-1)/4==c))?1:0, vec_view.value(i,k)[c] );
        //check gradients
        arma::rowvec gr = vec_view.grad(i,k).row(c);
        if ((i-1) / 4 == c)
        { // gradient of nonzero component
          switch ((i-1)%4)
          {
            case 0:
              EXPECT_ARMA_EQ( arma::rowvec("-1 -1 -1"), gr );
              break;
            case 1:
              EXPECT_ARMA_EQ( arma::rowvec("1 0 0"), gr );
              break;
            case 2:
              EXPECT_ARMA_EQ( arma::rowvec("0 1 0"), gr );
              break;
            case 3:
              EXPECT_ARMA_EQ( arma::rowvec("0 0 1"), gr );
              break;
            case 4: // last basis function does not contribute to the vector
              EXPECT_ARMA_EQ( arma::rowvec("0 0 0"), gr );
              break;
          }
        }
        else
          EXPECT_ARMA_EQ( arma::rowvec("0 0 0"), gr );
      }
    
    // check value and gradient of P0 function
    EXPECT_EQ( 1, scalar_view.value(0,k) );
    EXPECT_ARMA_EQ( arma::vec("0 0 0"), scalar_view.grad(0,k) );
    
    // check RT0 function
    unsigned int dof_offset = fe_vec.n_dofs() + 1;
    for (unsigned int i=0; i<fe_sys.n_dofs(); i++)
    {
      arma::vec exp_value;
      switch (i-dof_offset)
      {
        case 0:
          exp_value = q.point(k)*2 - arma::vec("0 0 2");
          break;
        case 1:
          exp_value = q.point(k)*2 - arma::vec("0 2 0");
          break;
        case 2:
          exp_value = q.point(k)*2 - arma::vec("2 0 0");
          break;
        case 3:
          exp_value = q.point(k)*2;
          break;
        default:
          exp_value = arma::vec("0 0 0");
      }
      EXPECT_ARMA_EQ( exp_value, rt_view.value(i,k) );
    }
  }
}



