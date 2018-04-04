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

#include "fem/xfe_values.hh"
#include "fem/singularity.hh"
#include "fem/xfem_element_data.hh"
#include "fem/fe_rt.hh"
#include "quadrature/qxfem_factory.hh"
#include "quadrature/qxfem.hh"


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
    FE_P_disc<dim> fe(0);
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


void create_xfem_data(ElementFullIter ele,
                      std::shared_ptr<Singularity<0>> sing,
                      XFEMElementSingularData<2>* &xdata)
{
    cout << "create\n" << endl;
    //create singularity
    Space<3>::Point center ({0.5, 0.5, 0});
    Space<3>::Point direction_vector ({0,0,1});
    Space<3>::Point n = arma::cross(ele->node[1]->point() - ele->node[0]->point(),
                                    ele->node[2]->point() - ele->node[0]->point());
    sing = std::make_shared<Singularity<0>>(center, 0.03, direction_vector, n, 1000);
    sing->set_sigma(10);
    sing->set_pressure(100);
    
    //create xfem data
    xdata = new XFEMElementSingularData<2>();
    //TODO: set number of quantities
    xdata->global_enriched_dofs().resize(2);
    xdata->global_enriched_dofs()[0].resize(1);
    xdata->global_enriched_dofs()[1].resize(1);
    xdata->set_element(0, 0); // 2d ele index, 1d ele index
    xdata->add_data(sing, 0);   // sing index
    ele->xfem_data = xdata;
    
    std::vector<std::vector<int>>& dofs = xdata->global_enriched_dofs()[1];
    cout << xdata->n_enrichments() << endl;
    cout << ele->n_nodes() << endl;
//     dofs.resize(xdata->n_enrichments(), std::vector<int>(ele->n_nodes(), -1));.
    dofs.resize(xdata->n_enrichments());
    for(unsigned int w=0; w < xdata->n_enrichments(); w++){
        dofs[w].resize(ele->n_nodes(), -1);
        cout << dofs[w].size() << endl;
        for(unsigned int i=0; i < ele->n_nodes(); i++){
                cout << w << "  " << i << endl;
                dofs[w][i] = w*i + i;
            }
        }
//     xdata.print(cout);
    ele->xfem_data = xdata;
}

TEST(FeValues, test_all) {
  // integrate a polynomial defined on the ref. element over an arbitrary element

    {
        // 2d case: triangle (0,1) (2,0) (3,4) surface = 3*4 - 1*2/2 - 1*4/4 - 3*3/2 = 9/2, det(jac) = 9
        NodeVector nodes(3);
        nodes.add_item(0);
        nodes[0].point()[0] = 0.0;
        nodes[0].point()[1] = 0.0;
        nodes[0].point()[2] = 0.0;

        nodes.add_item(1);
        nodes[1].point()[0] = 2.0;
        nodes[1].point()[1] = 0.0;
        nodes[1].point()[2] = 0.0;

        nodes.add_item(2);
        nodes[2].point()[0] = 0.0;
        nodes[2].point()[1] = 2.0;
        nodes[2].point()[2] = 0.0;

        ElementVector el_vec(1);
        el_vec.add_item(0);

        RegionIdx reg;
        Element ele(2, NULL, reg);      //NULL - mesh pointer, empty RegionIdx

        ele.node = new Node * [ele.n_nodes()];
        for(int i =0; i < 3; i++) ele.node[i] = nodes(i);
        el_vec[0] = ele; // dangerous since Element has no deep copy constructor.

        ElementFullIter it( el_vec(0) );
        EXPECT_DOUBLE_EQ( 6+1.0/3, integrate<2>( it ) );

        std::shared_ptr<Singularity<0>> sing;
        XFEMElementSingularData<2> *xdata = nullptr;
        create_xfem_data(it, sing, xdata);
        
        QXFEMFactory* qfact = new QXFEMFactory();
        std::shared_ptr<QXFEM<2,3>> qxfem = qfact->create_singular(xdata->sing_vec(), it);
        
        delete qfact;
//         // projection methods
//         MappingP1<2,3> mapping;
//         arma::mat::fixed<3, 3> map = mapping.element_map(ele);
//         EXPECT_ARMA_EQ( arma::mat("0 2 3; 1 0 4; 0 0 0"), map);
//         EXPECT_ARMA_EQ( arma::vec("0.6 0.2 0.2"), mapping.project_real_to_unit( arma::vec("1.0 1.4 0.0"), map ) );
        
        FE_P_disc<2> fe_p0(0);
        FE_RT0<2> fe;
        MappingP1<2,3> map;
        cout << "create xfe_values\n" << endl;
        XFEValues<2,3> xfe_values(map, fe, fe_p0, update_JxW_values | update_quadrature_points);
        
        xfe_values.reinit(it, *xdata, *qxfem);
        
        EXPECT_EQ(qxfem->size(), xfe_values.n_points());
        EXPECT_EQ(4, xfe_values.n_dofs());
        
        
        cout << "create xfe_values_p\n" << endl;
        XFEValues<2,3> xfe_values_p(map, fe_p0, fe_p0, update_JxW_values | update_quadrature_points);
        
        xfe_values_p.reinit(it, *xdata, *qxfem);
        
        EXPECT_EQ(qxfem->size(), xfe_values_p.n_points());
        EXPECT_EQ(2, xfe_values_p.n_dofs());
        
        if(xdata != nullptr)
            delete xdata;
    }

}
