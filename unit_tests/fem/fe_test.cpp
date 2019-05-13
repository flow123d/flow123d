#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest.hh>
#include "fem/fe_p.hh"                   // for FE_p, FE_CR
#include "fem/fe_rt.hh"                  // for FE_RT0
#include "mesh/ref_element.hh"           // for RefElement
#include "quadrature/quadrature_lib.hh"  // for QGauss



template<unsigned int dim>
void test_p()
{
    // test polynomial degree up to 5
    for (unsigned degree=0; degree<5; degree++)
    {
        FE_P<dim> fe_p(degree);
        
        // test number of components
        EXPECT_EQ( fe_p.n_components(), 1 );
        
        // test number of dofs
        switch (dim)
        {
            case 0:
                EXPECT_EQ( fe_p.n_dofs(), 1 );
                break;
            case 1:
                EXPECT_EQ( fe_p.n_dofs(), degree+1 );
                break;
            case 2:
                EXPECT_EQ( fe_p.n_dofs(), (degree+1)*(degree+2)/2 );
                break;
            case 3:
                EXPECT_EQ( fe_p.n_dofs(), (degree+1)*(degree*degree+5*degree+6)/6 );
                break;
        }
        
        // test partition of unity at quadrature points
        QGauss<dim> q(5);
        for (unsigned int k=0; k<q.size(); k++)
        {
            double sum = 0;
            for (unsigned int idof=0; idof<fe_p.n_dofs(); idof++)
                sum += fe_p.shape_value(idof, q.point(k), 0);
            
            EXPECT_NEAR( sum, 1, 2e-14 );
        }
    }
}


template<unsigned int dim>
void test_cr()
{
    FE_CR<dim> fe_cr;
    
    // test number of components
    EXPECT_EQ( fe_cr.n_components(), 1 );
    
    // test number of dofs
    EXPECT_EQ( fe_cr.n_dofs(), dim+1 );
    
    // test shape values at side barycenters
    arma::vec::fixed<dim> p;
    for (unsigned int sid=0; sid<dim+1; sid++)
    {
        // compute coordinates of side barycenter
        p.zeros();
        for (unsigned i=0; i<dim; i++)
            p += RefElement<dim>::node_coords(RefElement<dim>::interact(Interaction<0,dim-1>(sid))[i]);
        p /= dim;
        
        // shape values should equal 0 or 1
        for (unsigned int idof=0; idof<dim+1; idof++)
            EXPECT_DOUBLE_EQ( fe_cr.shape_value(idof, p, 0), idof==sid?1:0 );
    }
}


template<unsigned int dim>
void test_rt0()
{
    FE_RT0<dim> fe_rt;
    
    // test number of components
    EXPECT_EQ( fe_rt.n_components(), dim );
    
    // test number of dofs
    EXPECT_EQ( fe_rt.n_dofs(), dim+1 );
    
    // test fluxes at side barycenters
    arma::vec::fixed<dim> p;
    for (unsigned int sid=0; sid<dim+1; sid++)
    {
        // compute coordinates of side barycenter
        p.zeros();
        for (unsigned i=0; i<dim; i++)
            p += RefElement<dim>::node_coords(RefElement<dim>::interact(Interaction<0,dim-1>(sid))[i]);
        p /= dim;
        
        // fluxes should equal 0 o 1
        for (unsigned int idof=0; idof<dim+1; idof++)
        {
            double flux = 0;
            for (unsigned int c=0; c<dim; c++)
                flux += fe_rt.shape_value(idof, p, c) * RefElement<dim>::normal_vector(sid)[c];
            EXPECT_DOUBLE_EQ( flux*RefElement<dim>::side_measure(sid), idof==sid?1:0 );
        }
    }
}



TEST(FETest, FE_P) {
    test_p<0>();
    test_p<1>();
    test_p<2>();
    test_p<3>();
}


TEST(FETest, FE_CR) {
    test_cr<0>();
    test_cr<1>();
    test_cr<2>();
    test_cr<3>();
}


TEST(FETest, FE_RT0) {
//     test_rt0<0>();
    test_rt0<1>();
    test_rt0<2>();
    test_rt0<3>();
}



