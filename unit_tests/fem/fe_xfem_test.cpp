/*
 * 
 * 
 */

#define TEST_USE_PETSC
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>

#include "mesh/point.hh"
#include "mesh/elements.h"
#include "mesh/mesh.h"
#include "mesh_constructor.hh"

#include "fem/singularity.hh"
#include "fem/fe_rt.hh"
#include "fem/fe_rt_xfem.hh"
#include "fem/fe_p0_xfem.hh"
#include "fem/fe_values.hh"
#include "fem/mapping_p1.hh"

#include "quadrature/quadrature_lib.hh"
#include "quadrature/qxfem.hh"
#include "quadrature/qxfem_factory.hh"

typedef Space<3>::Point Point;

// simplest mesh
string ref_element_mesh = R"CODE(
$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
3
1 0 0 0
2 3 0 2
3 0 6 4
$EndNodes
$Elements
1
1 2 2 39 40 1 3 2
$EndElements
)CODE";

void print_fv_side(ElementFullIter ele, FESideValues<2,3>& fv_side, unsigned int side){
    fv_side.reinit(ele, side);
    
    cout << "FV_SIDE: side " << side << endl;
    cout << "n qpoints " << fv_side.n_points() << endl;
    cout << "ndofs " << fv_side.n_dofs() << endl;
//     cout << fv_side.shape_vector(0,0) << endl;
//     cout << fv_side.normal_vector(0) << endl;
//     cout << side_quad.weight(0) << endl;
    cout << "side measure " << ele->side(side)->measure() << endl;
    cout << arma::dot(fv_side.shape_vector(0,0),fv_side.normal_vector(0))*ele->side(side)->measure() << endl;
    cout << arma::dot(fv_side.shape_vector(1,0),fv_side.normal_vector(0))*ele->side(side)->measure() << endl;
    cout << arma::dot(fv_side.shape_vector(2,0),fv_side.normal_vector(0))*ele->side(side)->measure() << endl;
    cout << arma::dot(fv_side.shape_vector(3,0),fv_side.normal_vector(0))*ele->side(side)->measure() << endl;
    cout << arma::dot(fv_side.shape_vector(4,0),fv_side.normal_vector(0))*ele->side(side)->measure() << endl;
    cout << arma::dot(fv_side.shape_vector(5,0),fv_side.normal_vector(0))*ele->side(side)->measure() << endl;
    cout << endl;
    
    double sv_norm = arma::norm(fv_side.shape_vector(3,0),2);
    arma::vec sx = fv_side.shape_vector(3,0);
    if(sv_norm != 0) sx = sx / sv_norm;
    
    arma::vec n = fv_side.normal_vector(0) / arma::norm(fv_side.normal_vector(0),2);
    fv_side.shape_vector(3,0).print(cout,"shape");
    sx.print(cout,"sx");
    n.print(cout,"n");
    cout << "sx.n = " << arma::dot(sx,n) << endl;
    
    
}

TEST(fe_xfem, fe_rt_xfem) {

    QXFEMFactory<2,3> qfactory;
    
    // read mesh - simplset cube from test1
    Mesh* mesh = mesh_constructor();
    stringstream in(ref_element_mesh.c_str());
    mesh->read_gmsh_from_stream(in);
    ElementFullIter ele = mesh->element(0);
    Point n = arma::cross(ele->node[1]->point() - ele->node[0]->point(),
                          ele->node[2]->point() - ele->node[0]->point());
    
    auto func = std::make_shared<Singularity0D<3>>(arma::vec({1,2,2}),0.1,arma::vec({0,0,1}),n);
    shared_ptr<QXFEM<2,3>> qxfem = qfactory.create_singular({func},ele);
    
//     string dir_name = string(UNIT_TESTS_SRC_DIR) + "/fem/qxfem_output/";
//     qfactory.gnuplot_refinement(ele, dir_name, *qxfem2, {func});

    
//     std::ofstream q_points_file;
//     q_points_file.open (dir_name + "unit_q_points.dat");
//     if (q_points_file.is_open()) 
//     {
//         for(const Space<2>::Point &p : qxfem2->get_points())
//             q_points_file << p[0] << " " << p[1] << " " << p[2] << "\n";
//     }
//     else 
//     { 
//         std::cout << "Coud not write refinement for gnuplot.\n";
//     }
//     q_points_file.close();
        
//     func.evaluate_q_points(100, ele);
//     std::ofstream q_points_file;
//     q_points_file.open (dir_name + "q_points.dat");
//     if (q_points_file.is_open()) 
//     {
//         for(const Singularity0D<3>::Point &p : func.q_points())
//         q_points_file << p[0] << " " << p[1] << " " << p[2] << "\n";
//     }
//     else 
//     { 
//         std::cout << "Coud not write refinement for gnuplot.\n";
//     }
//     q_points_file.close();
    // add this to splot: 'q_points.dat' using 1:2:3 with points lc rgb 'green' title 'q_points'   
    
    MappingP1<2,3> map;
    
    FE_RT0<2,3> fe_rt;
    FE_RT0_XFEM<2,3> fe_rt_xfem(&fe_rt,{func});
    FEValues<2,3> fe_values(map, *qxfem, fe_rt_xfem, update_values | update_JxW_values | update_jacobians | update_inverse_jacobians | update_quadrature_points | update_divergence);
    
    fe_values.reinit(ele);
    cout << fe_values.n_points() << endl;
    cout << fe_values.n_dofs() << endl;
    cout << fe_values.determinant(0) << endl;
    cout << fe_values.shape_vector(3,0)[0] << "  " << fe_values.shape_vector(3,0)[1] << endl;
    cout << fe_values.shape_vector(5,0)[0] << "  " << fe_values.shape_vector(5,0)[1] << endl;
    cout << fe_values.shape_divergence(0,0) << endl;
    cout << fe_values.shape_divergence(1,0) << endl;
    cout << fe_values.shape_divergence(2,0) << endl;
    cout << fe_values.shape_divergence(3,0) << endl;
    cout << fe_values.shape_divergence(4,0) << endl;
    cout << fe_values.shape_divergence(5,0) << endl;
//     cout << fe_values.shape_vector(6,0)[0] << "  " << fe_values.shape_vector(6,0)[1] << endl;
    
    
    FE_P_disc<0,2,3> fe_p0;
    FE_P0_XFEM<2,3> fe_p0_xfem(&fe_p0,{func});
    FEValues<2,3> fe_values_p(map, *qxfem, fe_p0_xfem, update_values | update_JxW_values | update_inverse_jacobians | update_quadrature_points);
    
    fe_values_p.reinit(ele);
    
    cout << "Side values:" << endl;
    QGauss<1> side_quad(1);
    FESideValues<2,3> fv_side(map, side_quad, fe_rt_xfem, update_values | update_JxW_values
                                                          | update_quadrature_points | update_normal_vectors);
    
    print_fv_side(ele,fv_side, 0);
    print_fv_side(ele,fv_side, 1);
    print_fv_side(ele,fv_side, 2);
//     fv_side.reinit(ele, 0);
//     cout << fv_side.n_points() << endl;
//     cout << fv_side.n_dofs() << endl;
//     cout << fv_side.shape_vector(0,0) << endl;
//     cout << fv_side.normal_vector(0) << endl;
//     cout << side_quad.weight(0) << endl;
//     cout << ele->side(0)->measure() << endl;
//     cout << arma::dot(fv_side.shape_vector(0,0),fv_side.normal_vector(0))*ele->side(0)->measure() << endl;
//     cout << arma::dot(fv_side.shape_vector(1,0),fv_side.normal_vector(0))*ele->side(0)->measure() << endl;
//     cout << arma::dot(fv_side.shape_vector(2,0),fv_side.normal_vector(0))*ele->side(0)->measure() << endl;
//     cout << arma::dot(fv_side.shape_vector(3,0),fv_side.normal_vector(0))*ele->side(0)->measure() << endl;
//     cout << arma::dot(fv_side.shape_vector(4,0),fv_side.normal_vector(0))*ele->side(0)->measure() << endl;
//     cout << arma::dot(fv_side.shape_vector(5,0),fv_side.normal_vector(0))*ele->side(0)->measure() << endl;
    
}

