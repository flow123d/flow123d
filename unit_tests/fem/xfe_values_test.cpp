#define FEAL_OVERRIDE_ASSERTS

#include <flow_gtest.hh>
#include <cmath>
#include <armadillo>
#include "arma_expect.hh"

#include "system/armadillo_tools.hh"
#include "system/sys_profiler.hh"
#include "quadrature/quadrature_lib.hh"
#include "fem/fe_p.hh"
#include "fem/fe_values.hh"
#include "fem/mapping_p1.hh"
#include "mesh/mesh.h"
#include "mesh/accessors.hh"
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
double integrate(ElementAccessor<3> &ele) {
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


void create_xfem_data_2d(ElementAccessor<3> &ele,
                      std::shared_ptr<Singularity<0>> &sing,
                      XFEMElementSingularData<2>* &xdata)
{
    cout << "create 2d xfem data\n" << endl;
    //create singularity
    Space<3>::Point center ({0.5, 0.5, 0});
    Space<3>::Point direction_vector ({0,0,1});
    Space<3>::Point n = arma::cross(ele.node(1)->point() - ele.node(0)->point(),
                                    ele.node(2)->point() - ele.node(0)->point());
    sing = std::make_shared<Singularity<0>>(center, 0.03, direction_vector, n, 1000);
    sing->set_sigma(10);
    sing->set_pressure(100);
    
    //create xfem data
    xdata = new XFEMElementSingularData<2>();
    //TODO: set number of quantities
    xdata->global_enriched_dofs().resize(2);
    xdata->global_enriched_dofs()[0].resize(1);
    xdata->global_enriched_dofs()[1].resize(1);
    xdata->set_element(0); // 2d ele index
    xdata->add_data(sing, 0, 0);   // sing index, 1d ele index
    
    std::vector<std::vector<int>>& dofs = xdata->global_enriched_dofs()[1];
//     dofs.resize(xdata->n_enrichments(), std::vector<int>(ele->n_nodes(), -1));.
    dofs.resize(xdata->n_enrichments());
    for(unsigned int w=0; w < xdata->n_enrichments(); w++){
        dofs[w].resize(ele->n_nodes(), -1);
        for(unsigned int i=0; i < ele->n_nodes(); i++){
                dofs[w][i] = w*i + i;
            }
        }
//     xdata.print(cout);
}

TEST(FeValues, test_2d) {
    {
        Mesh mesh;
        mesh.init_node_vector(3);
        mesh.add_node(0, {0,0,0});
        mesh.add_node(1, {2,0,0});
        mesh.add_node(2, {0,2,0});
        
        // elm_id, int dim, region_id, partition_id, std::vector<unsigned int> node_ids
        mesh.init_element_vector(1);
        mesh.add_element(0, 2, 0, 0, {0,1,2});

        ElementAccessor<3> ele(&mesh, 0);

        mesh.side_nodes.resize(ele->n_nodes());
        mesh.side_nodes[1] = {{0,1},{0,2},{2,1}};

        EXPECT_DOUBLE_EQ( 6+1.0/3, integrate<2>( ele ) );

        std::shared_ptr<Singularity<0>> sing;
        XFEMElementSingularData<2> *xdata = nullptr;
        create_xfem_data_2d(ele, sing, xdata);
        
        QXFEMFactory qfact;
        std::shared_ptr<QXFEM<2,3>> qxfem = qfact.create_singular(xdata->sing_vec(), ele);
        
        FE_P_disc<2> fe_p0(0);
        FE_RT0<2> fe;
        MappingP1<2,3> map;
        XFEValues<2,3> xfe_values(map, fe, fe_p0, update_values |
                                                    update_JxW_values | update_jacobians |
                                                    update_inverse_jacobians | update_quadrature_points
                                                    | update_divergence);
        xfe_values.reinit(ele, *xdata, *qxfem);
        auto velocity = xfe_values.vector_view(0);
        
        EXPECT_EQ(qxfem->size(), xfe_values.n_points());
        EXPECT_EQ(4, xfe_values.n_dofs());
        {
            double sum = 0,
                   sum_enr = 0;
            for(unsigned int q=0; q < xfe_values.n_points(); q++) {
                sum += arma::dot(velocity.value(0,q),velocity.value(0,q)) * xfe_values.JxW(q);
                sum_enr += arma::dot(velocity.value(3,q),velocity.value(3,q)) * xfe_values.JxW(q);
            }
            EXPECT_NEAR(0.3328906523792483, sum, 1e-12);
            EXPECT_NEAR(459.220342181639, sum_enr, 1e-12);
        }
        
        XFEValues<2,3> xfe_values_p(map, fe_p0, fe_p0, update_values |
                                                        update_JxW_values | update_jacobians |
                                                        update_inverse_jacobians | update_quadrature_points
                                                        | update_divergence);
        xfe_values_p.reinit(ele, *xdata, *qxfem);
        auto pressure = xfe_values_p.scalar_view(0);
        
        EXPECT_EQ(qxfem->size(), xfe_values_p.n_points());
        EXPECT_EQ(2, xfe_values_p.n_dofs());
        
        {
            double sum = 0.0;
            for(unsigned int q=0; q < xfe_values_p.n_points(); q++) {
                sum += pressure.value(0,q) * xfe_values_p.JxW(q);
            }
//             cout << ele.measure() - sing->geometry().volume()- sum << endl;
            EXPECT_NEAR(ele.measure() - sing->geometry().volume(), sum, 6e-6);
        }
        
        if(xdata != nullptr)
            delete xdata;
    }
}



void create_xfem_data_3d(ElementAccessor<3> &ele,
                      std::shared_ptr<Singularity<1>> &sing,
                      XFEMElementSingularData<3>* &xdata)
{
    cout << "create 3d xfem data\n" << endl;
    //create singularity
    Space<3>::Point a ({0.5, 0.5, 0});
    Space<3>::Point b ({0.5, 0.5, 2});
    
    sing = std::make_shared<Singularity<1>>(a, b, 0.03, 1000, 100);
    sing->set_sigma(10);
    sing->set_pressure(100);
    
    //create xfem data
    xdata = new XFEMElementSingularData<3>();
    //TODO: set number of quantities
    xdata->global_enriched_dofs().resize(2);
    xdata->global_enriched_dofs()[0].resize(1);
    xdata->global_enriched_dofs()[1].resize(1);
    xdata->set_element(0); // 3d ele index, 1d ele index
    xdata->add_data(sing, 0, 0);   // sing index
    
    std::vector<std::vector<int>>& dofs = xdata->global_enriched_dofs()[1];
//     dofs.resize(xdata->n_enrichments(), std::vector<int>(ele->n_nodes(), -1));.
    dofs.resize(xdata->n_enrichments());
    for(unsigned int w=0; w < xdata->n_enrichments(); w++){
        dofs[w].resize(ele->n_nodes(), -1);
        for(unsigned int i=0; i < ele->n_nodes(); i++){
                dofs[w][i] = w*i + i;
            }
        }
//     xdata.print(cout);
}

TEST(FeValues, test_3d) {
        Mesh mesh;
        mesh.init_node_vector(4);
        mesh.add_node(0, {0,0,0});
        mesh.add_node(1, {2,0,0});
        mesh.add_node(2, {0,2,0});
        mesh.add_node(3, {0,0,2});
        
        // elm_id, int dim, region_id, partition_id, std::vector<unsigned int> node_ids
        mesh.init_element_vector(1);
        mesh.add_element(0, 3, 0, 0, {0,1,2,3});

        ElementAccessor<3> ele(&mesh, 0);

        mesh.side_nodes.resize(ele->n_nodes());
        mesh.side_nodes[2] = {{ 0, 1, 2 }, { 0, 1, 3 }, { 0, 2, 3 }, { 1, 2, 3 }};

        std::shared_ptr<Singularity<1>> sing;
        XFEMElementSingularData<3> *xdata = nullptr;
        create_xfem_data_3d(ele, sing, xdata);
        
        QXFEMFactory qfact(7);
        std::shared_ptr<QXFEM<3,3>> qxfem = qfact.create_singular(xdata->sing_vec(), ele);
        
        FE_P_disc<3> fe_p0(0);
        FE_RT0<3> fe;
        MappingP1<3,3> map;
        XFEValues<3,3> xfe_values(map, fe, fe_p0, update_values |
                                                    update_JxW_values | update_jacobians |
                                                    update_inverse_jacobians | update_quadrature_points
                                                    | update_divergence);
        xfe_values.reinit(ele, *xdata, *qxfem);
        auto velocity = xfe_values.vector_view(0);
        
        EXPECT_EQ(qxfem->size(), xfe_values.n_points());
        EXPECT_EQ(5, xfe_values.n_dofs());
        
        {
            double sum = 0,
                sum_enr = 0;
            for(unsigned int q=0; q < xfe_values.n_points(); q++) {
                sum += arma::dot(velocity.value(0,q),velocity.value(0,q)) * xfe_values.JxW(q);
                sum_enr += arma::dot(velocity.value(4,q),velocity.value(4,q)) * xfe_values.JxW(q);
            }
//             cout << setprecision(15) << sum << endl;
//             cout << setprecision(15) << sum_enr << endl;
            EXPECT_NEAR(0.26616240907934, sum, 1e-12);
            EXPECT_NEAR(119.995717102827, sum_enr, 1e-12);
        }
        
        XFEValues<3,3> xfe_values_p(map, fe_p0, fe_p0, update_values |
                                                        update_JxW_values | update_jacobians |
                                                        update_inverse_jacobians | update_quadrature_points
                                                        | update_divergence);
        xfe_values_p.reinit(ele, *xdata, *qxfem);
        auto pressure = xfe_values_p.scalar_view(0);
        
        EXPECT_EQ(qxfem->size(), xfe_values_p.n_points());
        EXPECT_EQ(2, xfe_values_p.n_dofs());
        
        {
            double sum = 0.0;
            for(unsigned int q=0; q < xfe_values_p.n_points(); q++) {
                sum += pressure.value(0,q) * xfe_values_p.JxW(q);
            }
//             cout << ele.measure() - sing->geometry().volume()/2 - sum << endl;
            EXPECT_NEAR(ele.measure() - sing->geometry().volume()/2, sum, 3e-5);
        }
        
        if(xdata != nullptr)
            delete xdata;
}
