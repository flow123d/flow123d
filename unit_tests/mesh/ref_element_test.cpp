/*
 * Ref Element unit test.
 * 
 * There is not much to test actually, but this serves mainly for user checking and testing.
 *
 * Author: pe
 */


#include <flow_gtest.hh>

#include "mesh/ref_element.hh"

#include "system/sys_profiler.hh"

using namespace std;

template<unsigned int dim>
void test_coordinates(){
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
void test_topology(){
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
    }
}

TEST(RefElement, test_coordinates) {
    Profiler::initialize();

    test_coordinates<1>();
    test_coordinates<2>();
    test_coordinates<3>();
    
    test_topology<1>();
    test_topology<2>();
    test_topology<3>();
}
