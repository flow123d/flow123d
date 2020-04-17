/*
 * Ref Element unit test.
 * 
 * ref_element_test.cpp
 *
 *  Created on: Jul 4, 2016
 *      Author: jb
 */



#include <flow_gtest.hh>
#include "mesh/ref_element.hh"

#include "../../src/system/armadillo_tools.hh"
#include "arma_expect.hh"
#include "system/sys_profiler.hh"
#include "armadillo"

using namespace std;

TEST(RefElement, barycentric_on_face) {
    armadillo_setup();

    // dim 1
    EXPECT_ARMA_EQ( arma::vec("1"),
            RefElement<1>::barycentric_on_face( arma::vec("1 0"), 0));
    EXPECT_ARMA_EQ( arma::vec("1"),
            RefElement<1>::barycentric_on_face( arma::vec("0 1"), 1));

    // dim 2
    EXPECT_ARMA_EQ( arma::vec("1 0"),             RefElement<2>::barycentric_on_face( arma::vec("1 0 0"), 0));
    EXPECT_ARMA_EQ( arma::vec("0 1"),             RefElement<2>::barycentric_on_face( arma::vec("0 1 0"), 0));

    EXPECT_ARMA_EQ( arma::vec("1 0"),             RefElement<2>::barycentric_on_face( arma::vec("1 0 0"), 1));
    EXPECT_ARMA_EQ( arma::vec("0 1"),             RefElement<2>::barycentric_on_face( arma::vec("0 0 1"), 1));

    EXPECT_ARMA_EQ( arma::vec("1 0"),             RefElement<2>::barycentric_on_face( arma::vec("0 1 0"), 2));
    EXPECT_ARMA_EQ( arma::vec("0 1"),             RefElement<2>::barycentric_on_face( arma::vec("0 0 1"), 2));


    // dim 3
    EXPECT_ARMA_EQ( arma::vec("1 0 0"),             RefElement<3>::barycentric_on_face( arma::vec("1 0 0 0"), 0));
    EXPECT_ARMA_EQ( arma::vec("0 1 0"),             RefElement<3>::barycentric_on_face( arma::vec("0 1 0 0"), 0));
    EXPECT_ARMA_EQ( arma::vec("0 0 1"),             RefElement<3>::barycentric_on_face( arma::vec("0 0 1 0"), 0));

    EXPECT_ARMA_EQ( arma::vec("1 0 0"),             RefElement<3>::barycentric_on_face( arma::vec("1 0 0 0"), 1));
    EXPECT_ARMA_EQ( arma::vec("0 1 0"),             RefElement<3>::barycentric_on_face( arma::vec("0 1 0 0"), 1));
    EXPECT_ARMA_EQ( arma::vec("0 0 1"),             RefElement<3>::barycentric_on_face( arma::vec("0 0 0 1"), 1));

    EXPECT_ARMA_EQ( arma::vec("1 0 0"),             RefElement<3>::barycentric_on_face( arma::vec("1 0 0 0"), 2));
    EXPECT_ARMA_EQ( arma::vec("0 1 0"),             RefElement<3>::barycentric_on_face( arma::vec("0 0 1 0"), 2));
    EXPECT_ARMA_EQ( arma::vec("0 0 1"),             RefElement<3>::barycentric_on_face( arma::vec("0 0 0 1"), 2));

    EXPECT_ARMA_EQ( arma::vec("1 0 0"),             RefElement<3>::barycentric_on_face( arma::vec("0 1 0 0"), 3));
    EXPECT_ARMA_EQ( arma::vec("0 1 0"),             RefElement<3>::barycentric_on_face( arma::vec("0 0 1 0"), 3));
    EXPECT_ARMA_EQ( arma::vec("0 0 1"),             RefElement<3>::barycentric_on_face( arma::vec("0 0 0 1"), 3));
}


TEST(RefElement, centers_of_subelements) {
    armadillo_setup();

    // dim = 1
    {
    auto list = RefElement<1>::centers_of_subelements(0);
    EXPECT_EQ( 2, list.size());
    EXPECT_ARMA_EQ( arma::vec("0"), list[0]);
    EXPECT_ARMA_EQ( arma::vec("1"), list[1]);
    }
    {
    auto list = RefElement<1>::centers_of_subelements(1);
    EXPECT_EQ( 1, list.size());
    EXPECT_ARMA_EQ( arma::vec("0.5"), list[0]);
    }

    // dim = 2
    {
    auto list = RefElement<2>::centers_of_subelements(0);
    EXPECT_EQ( 3, list.size());
    EXPECT_ARMA_EQ( arma::vec("0 0"), list[0]);
    EXPECT_ARMA_EQ( arma::vec("1 0"), list[1]);
    EXPECT_ARMA_EQ( arma::vec("0 1"), list[2]);
    }
    {
    auto list = RefElement<2>::centers_of_subelements(1);
    EXPECT_EQ( 3, list.size());
    EXPECT_ARMA_EQ( arma::vec("0.5 0"), list[0]);
    EXPECT_ARMA_EQ( arma::vec("0 0.5"), list[1]);
    EXPECT_ARMA_EQ( arma::vec("0.5 0.5"), list[2]);
    }
    {
    auto list = RefElement<2>::centers_of_subelements(2);
    EXPECT_EQ( 1, list.size());
    EXPECT_ARMA_EQ( arma::vec({ 1/3.0, 1/3.0 }), list[0]);
    }

}


TEST(RefElement, clip_1d) {
    armadillo_setup();

    // in element
    EXPECT_ARMA_EQ( arma::vec("0 1"), RefElement<1>::clip( arma::vec("0 1")));
    EXPECT_ARMA_EQ( arma::vec("0.5 0.5"), RefElement<1>::clip( arma::vec("0.5 0.5")));
    EXPECT_ARMA_EQ( arma::vec("1 0"), RefElement<1>::clip( arma::vec("1 0")));
    // out of element
    EXPECT_ARMA_EQ( arma::vec("0 1"), RefElement<1>::clip( arma::vec("-0.5 1.5")));
    EXPECT_ARMA_EQ( arma::vec("1 0"), RefElement<1>::clip( arma::vec("1.5 -0.5")));
}


TEST(RefElement, clip_2d) {
    armadillo_setup();

    //in element
    EXPECT_ARMA_EQ( arma::vec("1 0 0"), RefElement<2>::clip( arma::vec("1 0 0")));
    EXPECT_ARMA_EQ( arma::vec("0 1 0"), RefElement<2>::clip( arma::vec("0 1 0")));
    EXPECT_ARMA_EQ( arma::vec("0 0 1"), RefElement<2>::clip( arma::vec("0 0 1")));
    EXPECT_ARMA_EQ( arma::vec("0.3 0.3 0.4"), RefElement<2>::clip( arma::vec("0.3 0.3 0.4")));
    // out of element
    EXPECT_ARMA_EQ( arma::vec("0 0 1"), RefElement<2>::clip( arma::vec("-0.5 0 1.5")));
    EXPECT_ARMA_EQ( arma::vec("0 0 1"), RefElement<2>::clip( arma::vec("0 -0.5 1.5")));
    EXPECT_ARMA_EQ( arma::vec("0 1 0"), RefElement<2>::clip( arma::vec("-0.5 1 0.5")));
    EXPECT_ARMA_EQ( arma::vec("0 1 0"), RefElement<2>::clip( arma::vec("0.3 1.3 -0.6")));
    EXPECT_ARMA_EQ( arma::vec("1 0 0"), RefElement<2>::clip( arma::vec("1 -0.5 0.5")));
    EXPECT_ARMA_EQ( arma::vec("1 0 0"), RefElement<2>::clip( arma::vec("1.3 0.3 -0.6")));

    EXPECT_ARMA_EQ( arma::vec("0.5 0 0.5"), RefElement<2>::clip( arma::vec("0.5 -0.3 0.8")));
    EXPECT_ARMA_EQ( arma::vec("0 0.5 0.5"), RefElement<2>::clip( arma::vec("-0.3 0.5 0.8")));
    EXPECT_ARMA_EQ( arma::vec("0.5 0.5 0"), RefElement<2>::clip( arma::vec("1 1 -1")));
}


TEST(RefElement, interpolate) {
    armadillo_setup();
    
//     RefElement<1>::bary_coords<0>(0).print(cout,"1-0: 0");
//     RefElement<1>::bary_coords<0>(1).print(cout,"1-0: 1");
//     
//     RefElement<2>::bary_coords<1>(0).print(cout,"2-1: 0");
//     RefElement<2>::bary_coords<1>(1).print(cout,"2-1: 1");
//     RefElement<2>::bary_coords<1>(2).print(cout,"2-1: 2");
//     
//     RefElement<3>::bary_coords<1>(0).print(cout,"3-1: 0");
//     RefElement<3>::bary_coords<1>(1).print(cout,"3-1: 1");
//     RefElement<3>::bary_coords<1>(2).print(cout,"3-1: 2");
//     RefElement<3>::bary_coords<1>(3).print(cout,"3-1: 3");
//     RefElement<3>::bary_coords<1>(4).print(cout,"3-1: 4");
//     RefElement<3>::bary_coords<1>(5).print(cout,"3-1: 5");
//     
//     RefElement<3>::bary_coords<2>(0).print(cout,"3-2: 0");
//     RefElement<3>::bary_coords<2>(1).print(cout,"3-2: 1");
//     RefElement<3>::bary_coords<2>(2).print(cout,"3-2: 2");
//     RefElement<3>::bary_coords<2>(3).print(cout,"3-2: 3");
    
    //VF
//     EXPECT_ARMA_EQ( arma::vec("0.75 0.25 0"),       RefElement<2>::interpolate<1>("0.75 0.25",0));
//     EXPECT_ARMA_EQ( arma::vec("0.75 0 0.25"),       RefElement<2>::interpolate<1>("0.75 0.25",1));
//     EXPECT_ARMA_EQ( arma::vec("0 0.75 0.25"),       RefElement<2>::interpolate<1>("0.75 0.25",2));
//     
//     EXPECT_ARMA_EQ( arma::vec("0.75 0.25 0 0"),     RefElement<3>::interpolate<1>("0.75 0.25",0));
//     EXPECT_ARMA_EQ( arma::vec("0.75 0 0.25 0"),     RefElement<3>::interpolate<1>("0.75 0.25",1));
//     EXPECT_ARMA_EQ( arma::vec("0 0.75 0.25 0"),     RefElement<3>::interpolate<1>("0.75 0.25",2));
//     
//     EXPECT_ARMA_EQ( arma::vec("0.75 0 0 0.25"),     RefElement<3>::interpolate<1>("0.75 0.25",3));
//     EXPECT_ARMA_EQ( arma::vec("0 0.75 0 0.25"),     RefElement<3>::interpolate<1>("0.75 0.25",4));
//     EXPECT_ARMA_EQ( arma::vec("0 0 0.75 0.25"),     RefElement<3>::interpolate<1>("0.75 0.25",5));

    // dim 1
    EXPECT_ARMA_EQ( arma::vec("1 0"),               RefElement<1>::interpolate<0>( arma::vec("1"), 0));
    EXPECT_ARMA_EQ( arma::vec("0 1"),               RefElement<1>::interpolate<0>( arma::vec("1"), 1));

    // dim 2
    EXPECT_ARMA_EQ( arma::vec("1 0 0"),             RefElement<2>::interpolate<1>( arma::vec("1 0"), 0));
    EXPECT_ARMA_EQ( arma::vec("0 1 0"),             RefElement<2>::interpolate<1>( arma::vec("0 1"), 0));
    EXPECT_ARMA_EQ( arma::vec("1 0 0"),             RefElement<2>::interpolate<1>( arma::vec("1 0"), 1));
    EXPECT_ARMA_EQ( arma::vec("0 0 1"),             RefElement<2>::interpolate<1>( arma::vec("0 1"), 1));
    EXPECT_ARMA_EQ( arma::vec("0 1 0"),             RefElement<2>::interpolate<1>( arma::vec("1 0"), 2));
    EXPECT_ARMA_EQ( arma::vec("0 0 1"),             RefElement<2>::interpolate<1>( arma::vec("0 1"), 2));

    EXPECT_ARMA_EQ( arma::vec("0.25 0.75 0"),       RefElement<2>::interpolate<1>("0.25 0.75",0));
    EXPECT_ARMA_EQ( arma::vec("0.25 0 0.75"),       RefElement<2>::interpolate<1>("0.25 0.75",1));
    EXPECT_ARMA_EQ( arma::vec("0 0.25 0.75"),       RefElement<2>::interpolate<1>("0.25 0.75",2));
    
    // dim 3
    EXPECT_ARMA_EQ( arma::vec("1 0 0 0"),           RefElement<3>::interpolate<2>( arma::vec("1 0 0"), 0));
    EXPECT_ARMA_EQ( arma::vec("0 1 0 0"),           RefElement<3>::interpolate<2>( arma::vec("0 1 0"), 0));
    EXPECT_ARMA_EQ( arma::vec("0 0 1 0"),           RefElement<3>::interpolate<2>( arma::vec("0 0 1"), 0));

    EXPECT_ARMA_EQ( arma::vec("1 0 0 0"),           RefElement<3>::interpolate<2>( arma::vec("1 0 0"), 1));
    EXPECT_ARMA_EQ( arma::vec("0 1 0 0"),           RefElement<3>::interpolate<2>( arma::vec("0 1 0"), 1));
    EXPECT_ARMA_EQ( arma::vec("0 0 0 1"),           RefElement<3>::interpolate<2>( arma::vec("0 0 1"), 1));

    EXPECT_ARMA_EQ( arma::vec("1 0 0 0"),           RefElement<3>::interpolate<2>( arma::vec("1 0 0"), 2));
    EXPECT_ARMA_EQ( arma::vec("0 0 1 0"),           RefElement<3>::interpolate<2>( arma::vec("0 1 0"), 2));
    EXPECT_ARMA_EQ( arma::vec("0 0 0 1"),           RefElement<3>::interpolate<2>( arma::vec("0 0 1"), 2));

    EXPECT_ARMA_EQ( arma::vec("0 1 0 0"),           RefElement<3>::interpolate<2>( arma::vec("1 0 0"), 3));
    EXPECT_ARMA_EQ( arma::vec("0 0 1 0"),           RefElement<3>::interpolate<2>( arma::vec("0 1 0"), 3));
    EXPECT_ARMA_EQ( arma::vec("0 0 0 1"),           RefElement<3>::interpolate<2>( arma::vec("0 0 1"), 3));

    EXPECT_ARMA_EQ( arma::vec("0.25 0.75 0 0"),     RefElement<3>::interpolate<1>("0.25 0.75",0));
    EXPECT_ARMA_EQ( arma::vec("0.25 0 0.75 0"),     RefElement<3>::interpolate<1>("0.25 0.75",1));
    EXPECT_ARMA_EQ( arma::vec("0.25 0 0 0.75"),     RefElement<3>::interpolate<1>("0.25 0.75",2));
    EXPECT_ARMA_EQ( arma::vec("0 0.25 0.75 0"),     RefElement<3>::interpolate<1>("0.25 0.75",3));
    EXPECT_ARMA_EQ( arma::vec("0 0.25 0 0.75"),     RefElement<3>::interpolate<1>("0.25 0.75",4));
    EXPECT_ARMA_EQ( arma::vec("0 0 0.25 0.75"),     RefElement<3>::interpolate<1>("0.25 0.75",5));

    EXPECT_ARMA_EQ( arma::vec("0.5 0.3 0.2 0"),     RefElement<3>::interpolate<2>("0.5 0.3 0.2",0));
    EXPECT_ARMA_EQ( arma::vec("0.5 0.3 0 0.2"),     RefElement<3>::interpolate<2>("0.5 0.3 0.2",1));
    EXPECT_ARMA_EQ( arma::vec("0.5 0 0.3 0.2"),     RefElement<3>::interpolate<2>("0.5 0.3 0.2",2));
    EXPECT_ARMA_EQ( arma::vec("0 0.5 0.3 0.2"),     RefElement<3>::interpolate<2>("0.5 0.3 0.2",3));
    
}


TEST(RefElement, bary_local){
    arma::vec lp("0.3"), bp("0.7 0.3 ");
    EXPECT_ARMA_EQ( arma::vec(bp),  RefElement<1>::local_to_bary(lp));
    EXPECT_ARMA_EQ( arma::vec(lp),  RefElement<1>::bary_to_local(bp));
    
    lp = "0.2 0.3"; bp = "0.5 0.2 0.3";
    EXPECT_ARMA_EQ( arma::vec(bp),  RefElement<2>::local_to_bary(lp));
    EXPECT_ARMA_EQ( arma::vec(lp),  RefElement<2>::bary_to_local(bp));
    
    lp = "0.2 0.3 0.4"; bp = "0.1 0.2 0.3 0.4";
    EXPECT_ARMA_EQ( arma::vec(bp),  RefElement<3>::local_to_bary(lp));
    EXPECT_ARMA_EQ( arma::vec(lp),  RefElement<3>::bary_to_local(bp));
}








/*

// Write down the definition of RefElement

template<unsigned int dim>
void coordinates(){
    cout << "\ndim = " << dim  << "  ------------------\n" << endl;
    for(unsigned int nid=0; nid < RefElement<dim>::n_nodes; nid++)
    {
        arma::vec::fixed<dim> node = RefElement<dim>::node_coords(nid);
        cout << "node (local coordinates) " << nid << ":\t";
        for(unsigned int d=0; d < dim; d++)
            cout << node[d] << " ";
        cout << endl;
    }
    
    for(unsigned int nid=0; nid < RefElement<dim>::n_nodes; nid++)
    {
        arma::vec::fixed<dim+1> bnode = RefElement<dim>::node_barycentric_coords(nid);
        cout << "barycentric node " << nid << ":\t";
        for(unsigned int d=0; d < dim+1; d++)
            cout << bnode[d] << " ";
        cout << endl;
    }
    
    for(unsigned int sid=0; sid < RefElement<dim>::n_sides; sid++)
    {
        arma::vec::fixed<dim> normal = RefElement<dim>::normal_vector(sid);
        cout << "normal vector of side " << sid << ":\t";
        for(unsigned int d=0; d < dim; d++)
            cout << normal[d] << " ";
        cout << endl;
    }
    
    for(unsigned int sid=0; sid < RefElement<dim>::n_sides; sid++)
    {
        cout << "oposite node to side " << sid << ":\t" << RefElement<dim>::oposite_node(sid) 
        << "\t side measure: \t" << RefElement<dim>::side_measure(sid) << endl;
    }
    
    
    if(dim == 3)
    cout << "Jacobian = "
    << arma::dot( arma::cross(RefElement<3>::node_coords(1) - RefElement<3>::node_coords(0), 
                                  RefElement<3>::node_coords(2) - RefElement<3>::node_coords(0)),
                    RefElement<3>::node_coords(3) - RefElement<3>::node_coords(0)
                    )
    << endl;
}

template<unsigned int dim>
void topology(){
    cout << "\ndim = " << dim  << "  ------------------ topology\n" << endl;
    IdxVector<2> v({0,1});
    
    cout << "line nodes:\n";
    for(unsigned int i=0; i < RefElement<dim>::n_lines; i++)
    {
        for(unsigned int j=0; j < RefElement<1>::n_nodes; j++)
            cout << RefElement<dim>::template interact<0,1>(i)[j] << " ";
        cout << endl;
    }
    
    cout << "node lines:\n";
    for(unsigned int i=0; i < RefElement<dim>::n_nodes; i++)
    {
        for(unsigned int j=0; j < RefElement<dim>::n_lines_per_node; j++)
            cout << RefElement<dim>::template interact<1,0>(i)[j] << " ";
        cout << endl;
    }
    
    if(dim == 2)
    {
        cout << "topology_idx: dim=2 subdim=0 " << RefElement<2>::topology_idx<0>(3) << endl;
        cout << "topology_idx: dim=2 subdim=0 " << RefElement<2>::topology_idx<0>(5) << endl;
        cout << "topology_idx: dim=2 subdim=0 " << RefElement<2>::topology_idx<0>(6) << endl;
        cout << "topology_idx: dim=2 subdim=1 " << RefElement<2>::topology_idx<1>(1) << endl;
        cout << "topology_idx: dim=2 subdim=1 " << RefElement<2>::topology_idx<1>(2) << endl;
        cout << "topology_idx: dim=2 subdim=1 " << RefElement<2>::topology_idx<1>(4) << endl;
    }
    
    if(dim == 3)
    {
        cout << "side nodes:\n";
        for(unsigned int i=0; i < RefElement<dim>::n_sides; i++)
        {
            for(unsigned int j=0; j < RefElement<dim>::n_nodes_per_side; j++)
                cout << RefElement<dim>::template interact<0,2>(i)[j] << " ";
            cout << endl;
        }
        cout << "node sides:\n";
        for(unsigned int i=0; i < RefElement<dim>::n_nodes; i++)
        {
            for(unsigned int j=0; j < RefElement<dim>::n_nodes_per_side; j++)
                cout << RefElement<dim>::template interact<2,0>(i)[j] << " ";
            cout << endl;
        }
        
        cout << "line sides:\n";
        for(unsigned int i=0; i < RefElement<dim>::n_lines; i++)
        {
            for(unsigned int j=0; j < 2; j++)
                cout << RefElement<dim>::template interact<2,1>(i)[j] << " ";
            cout << endl;
        }
        
        cout << "side lines:\n";
        for(unsigned int i=0; i < RefElement<dim>::n_sides; i++)
        {
            for(unsigned int j=0; j < RefElement<dim>::n_lines_per_side; j++)
                cout << RefElement<dim>::template interact<1,2>(i)[j] << " ";
            cout << endl;
        }
        
        cout << "topology_idx: dim=3 subdim=2 " << RefElement<3>::topology_idx<2>(1) << endl;
        cout << "topology_idx: dim=3 subdim=2 " << RefElement<3>::topology_idx<2>(2) << endl;
        cout << "topology_idx: dim=3 subdim=2 " << RefElement<3>::topology_idx<2>(4) << endl;
        cout << "topology_idx: dim=3 subdim=2 " << RefElement<3>::topology_idx<2>(8) << endl;
        
        cout << "topology_idx: dim=3 subdim=1 " << RefElement<3>::topology_idx<1>(3) << endl;
        cout << "topology_idx: dim=3 subdim=1 " << RefElement<3>::topology_idx<1>(5) << endl;
        cout << "topology_idx: dim=3 subdim=1 " << RefElement<3>::topology_idx<1>(9) << endl;
        cout << "topology_idx: dim=3 subdim=1 " << RefElement<3>::topology_idx<1>(6) << endl;
        cout << "topology_idx: dim=3 subdim=1 " << RefElement<3>::topology_idx<1>(10) << endl;
        cout << "topology_idx: dim=3 subdim=1 " << RefElement<3>::topology_idx<1>(12) << endl;
        
        cout << "topology_idx: dim=3 subdim=0 " << RefElement<3>::topology_idx<0>(7) << endl;
        cout << "topology_idx: dim=3 subdim=0 " << RefElement<3>::topology_idx<0>(11) << endl;
        cout << "topology_idx: dim=3 subdim=0 " << RefElement<3>::topology_idx<0>(13) << endl;
        cout << "topology_idx: dim=3 subdim=0 " << RefElement<3>::topology_idx<0>(14) << endl;
    }
}

TEST(RefElement, write_down) {
    Profiler::instance();

    coordinates<1>();
    coordinates<2>();
    coordinates<3>();
    
    topology<1>();
    topology<2>();
    topology<3>();
}
*/
