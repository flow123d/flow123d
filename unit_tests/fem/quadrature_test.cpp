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

#define INTEGRATE(dim,  _func_ ) for( unsigned int i=0; i < quad.size(); i++) sum +=  _func_( quad.point<dim>(i)) * quad.weight(i);

double test_1_1d( const arma::vec::fixed<1> & p) {
    return 3 * p[0] + 1.0;
}

double test_2_1d( const arma::vec::fixed<1> & p) {
    return 3 * p[0] * p[0] + p[0] + 1.0;
}


TEST(Quadrature, test_1d) {

    {
    QGauss quad( 1, 1 ); // should integrate P1 exactly
    EXPECT_EQ(1, quad.size());
    double sum =0.0;
    INTEGRATE(1, test_1_1d);
    EXPECT_DOUBLE_EQ(5.0/2.0, sum); // 3 * 1/2 + 1
    }

    {
    QGauss quad( 1, 2 ); // should integrate P2 exactly
    EXPECT_EQ(2, quad.size());
    double sum =0.0;
    INTEGRATE(1, test_2_1d);
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
    QGauss quad( 2, 1 ); // should integrate P1 exactly
    EXPECT_EQ(1, quad.size());
    double sum =0.0;
    INTEGRATE(2,test_1_2d);
    EXPECT_DOUBLE_EQ(8.0/6.0, sum); // 3 * 1/6 + 2 * 1/6 + 1/2 = 8/6
    }

    {
    QGauss quad( 2, 2 ); // should integrate P2 exactly
    EXPECT_EQ(3, quad.size());
    double sum =0.0;
    INTEGRATE(2, test_2_2d);
    EXPECT_DOUBLE_EQ(19.0 / 12.0 , sum); // 3 * 1/12 + 1/6 + 6 * 1/12 + 1/6 + 1/2 = 19/12
    }
}


TEST(Quadrature, midpoint){
    QMidpoint quad(25);
    EXPECT_EQ(25, quad.size());
    
    double sum =0.0;
    INTEGRATE(1, test_1_1d);   // should integrate P1 exactly
    EXPECT_DOUBLE_EQ(5.0/2.0, sum); // 3 * 1/2 + 1
}


/// Map lower dimensional quadrature to element quadrature
/// using all sides and their permutations. Test equality of node coordinates
/// for both quadratures.
template<unsigned int dim>
void test_side_projection(Quadrature &subq)
{
	ASSERT_EQ(subq.dim(), dim-1);
    for (unsigned int sid=0; sid<RefElement<dim>::n_sides; sid++)
    {
        for (unsigned int pid=0; pid<RefElement<dim>::n_side_permutations; pid++)
        {
            Quadrature q = subq.make_from_side<dim>(sid, pid);
            
            std::vector<arma::vec::fixed<dim>> bary_subq;
            std::vector<arma::vec::fixed<dim+1>> bary_q;
            arma::mat::fixed<3,dim> coords_subq;
            arma::mat::fixed<3,dim+1> coords_q;
            
            // Setup coordinate matrices for mapping of reference points from side to 3d space.
            // We use the reference simplex as the 'spacial' element.
            for (unsigned int i=0; i<RefElement<dim>::n_nodes_per_side; i++)
            {
                unsigned int side_node_idx = RefElement<dim>::side_permutations[pid][i]; // index of permuted node within side
                unsigned int node_idx = RefElement<dim>::nodes_of_subelements[dim-1][sid][side_node_idx]; // index of node within element
                coords_subq.col(i) = arma::eye(3,dim)*RefElement<dim>::node_coords(node_idx);
            }
            // The same as above for mapping of reference element to 3d space.
            for (unsigned int i=0; i<RefElement<dim>::n_nodes; i++)
                coords_q.col(i) = arma::eye(3,dim)*RefElement<dim>::node_coords(i);

            // Setup barycentric coordinates of quadrature points.
            for (unsigned int i=0; i<subq.size(); i++)
            {
                bary_subq.push_back(RefElement<dim-1>::local_to_bary(subq.point<dim-1>(i)));
                bary_q.push_back(RefElement<dim>::local_to_bary(q.point<dim>(i)));
            }
            
            // Map barycentric coordinates of subquadrature and quadrature points to 3d space
            // and test their equality.
            for (unsigned int i=0; i<subq.size(); i++)
            {
                arma::vec3 sub_pts = coords_subq*bary_subq[i];
                arma::vec3 pts = coords_q*bary_q[i];
                EXPECT_ARMA_EQ( sub_pts, pts );
            }
        }
    }
}


TEST(Quadrature, side_projections){
    {
        for (unsigned int n_points=1; n_points<=10; n_points++)
        {
            QMidpoint qm(n_points);
            test_side_projection<2>(qm);
        }
    }
    
    {
        for (unsigned order=1; order<=10; order++)
        {
            QGauss qg(2, order);
            test_side_projection<3>(qg);
        }
    }
}

