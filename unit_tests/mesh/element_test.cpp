/*
 * element_test.cpp
 *
 *  Created on: Jul 6, 2016
 *      Author: jb
 */

#include <flow_gtest.hh>

#include "system/armadillo_setup.hh"
#include "arma_expect.hh"
#include "system/sys_profiler.hh"
#include "armadillo"

#include "mesh/element_impls.hh"
#include "mesh/region.hh"

class TestElement : public Element {
public:
    TestElement(std::vector<string> nodes_str)
    : Element()
    {
        std::vector<arma::vec3> nodes;
        for(auto str : nodes_str) nodes.push_back( arma::vec3(str));
        init(nodes.size()-1, nullptr, RegionIdx());
        unsigned int i=0;
        for(auto node : nodes)
            this->node[i++] = new Node(node[0], node[1], node[2]);
    }
};

TEST(Element, element_map) {
    Profiler::initialize();
    armadillo_setup();

    {
    TestElement ele({ "0 0 0", "1 0 0", "0 1 0", "0 0 1"});
    EXPECT_ARMA_EQ( arma::mat("1 0 0 0; 0 1 0 0; 0 0 1 0"), ele.element_map());
    EXPECT_ARMA_EQ( arma::vec("0.1 0.2 0.3 0.4"), ele.project_point( arma::vec3("0.1 0.2 0.3") ) );
    EXPECT_ARMA_EQ( arma::vec("0.5 0.5 0.5 -0.5"), ele.project_point( arma::vec3("0.5 0.5 0.5") ) );
    }

    {
    // trnaslated
    TestElement ele({ "1 2 3", "2 2 3", "1 3 3", "1 2 4"});
    EXPECT_ARMA_EQ( arma::mat("1 0 0 1; 0 1 0 2; 0 0 1 3"), ele.element_map());
    EXPECT_ARMA_EQ( arma::vec("0.1 0.2 0.3 0.4"), ele.project_point( arma::vec3("1.1 2.2 3.3") ) );
    EXPECT_ARMA_EQ( arma::vec("0.5 0.5 0.5 -0.5"), ele.project_point( arma::vec3("1.5 2.5 3.5") ) );
    }

    {
    // simplest cube element 7
    TestElement ele({ "-1 -1 1", "1 1 -1", "-1 -1 -1", "1 -1 -1"});
    EXPECT_ARMA_EQ( arma::mat("2 0 2 -1; 2 0 0 -1; -2 -2 -2 1"), ele.element_map());
    EXPECT_ARMA_EQ( arma::vec("0.25 0.25 0.25 0.25"), ele.project_point( arma::vec3("0 -0.5 -0.5") ) );
    //EXPECT_ARMA_EQ( arma::vec("0.1 0.2 0.3 0.4"), ele.project_point( arma::vec3("0.1 0.2 0.3") ) );
    }
}
