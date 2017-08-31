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
//     cout << arma::dot(fv_side.shape_vector(4,0),fv_side.normal_vector(0))*ele->side(side)->measure() << endl;
//     cout << arma::dot(fv_side.shape_vector(5,0),fv_side.normal_vector(0))*ele->side(side)->measure() << endl;
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

// TEST(fe_xfem, fe_rt_xfem) {
// 
//     QXFEMFactory<2,3> qfactory;
//     
//     // read mesh - simplset cube from test1
//     Mesh* mesh = mesh_constructor();
//     stringstream in(ref_element_mesh.c_str());
//     mesh->read_gmsh_from_stream(in);
//     ElementFullIter ele = mesh->element(0);
//     Point n = arma::cross(ele->node[1]->point() - ele->node[0]->point(),
//                           ele->node[2]->point() - ele->node[0]->point());
//     
//     auto func = std::make_shared<Singularity0D<3>>(arma::vec({1,2,2}),0.1,arma::vec({0,0,1}),n);
//     shared_ptr<QXFEM<2,3>> qxfem = qfactory.create_singular({func},ele);
//     
// //     string dir_name = string(UNIT_TESTS_SRC_DIR) + "/fem/qxfem_output/";
// //     qfactory.gnuplot_refinement(ele, dir_name, *qxfem2, {func});
// 
//     
// //     std::ofstream q_points_file;
// //     q_points_file.open (dir_name + "unit_q_points.dat");
// //     if (q_points_file.is_open()) 
// //     {
// //         for(const Space<2>::Point &p : qxfem2->get_points())
// //             q_points_file << p[0] << " " << p[1] << " " << p[2] << "\n";
// //     }
// //     else 
// //     { 
// //         std::cout << "Coud not write refinement for gnuplot.\n";
// //     }
// //     q_points_file.close();
//         
// //     func.evaluate_q_points(100, ele);
// //     std::ofstream q_points_file;
// //     q_points_file.open (dir_name + "q_points.dat");
// //     if (q_points_file.is_open()) 
// //     {
// //         for(const Singularity0D<3>::Point &p : func.q_points())
// //         q_points_file << p[0] << " " << p[1] << " " << p[2] << "\n";
// //     }
// //     else 
// //     { 
// //         std::cout << "Coud not write refinement for gnuplot.\n";
// //     }
// //     q_points_file.close();
//     // add this to splot: 'q_points.dat' using 1:2:3 with points lc rgb 'green' title 'q_points'   
//     
//     MappingP1<2,3> map;
//     
//     FE_RT0<2,3> fe_rt;
//     FE_RT0_XFEM<2,3> fe_rt_xfem(&fe_rt,{func});
//     FEValues<2,3> fe_values(map, *qxfem, fe_rt_xfem, update_values | update_JxW_values | update_jacobians | update_inverse_jacobians | update_quadrature_points | update_divergence);
//     
//     fe_values.reinit(ele);
//     cout << fe_values.n_points() << endl;
//     cout << fe_values.n_dofs() << endl;
//     cout << fe_values.determinant(0) << endl;
//     cout << fe_values.shape_vector(3,0)[0] << "  " << fe_values.shape_vector(3,0)[1] << endl;
//     cout << fe_values.shape_vector(5,0)[0] << "  " << fe_values.shape_vector(5,0)[1] << endl;
//     cout << fe_values.shape_divergence(0,0) << endl;
//     cout << fe_values.shape_divergence(1,0) << endl;
//     cout << fe_values.shape_divergence(2,0) << endl;
//     cout << fe_values.shape_divergence(3,0) << endl;
//     cout << fe_values.shape_divergence(4,0) << endl;
//     cout << fe_values.shape_divergence(5,0) << endl;
// //     cout << fe_values.shape_vector(6,0)[0] << "  " << fe_values.shape_vector(6,0)[1] << endl;
//     
//     
//     FE_P_disc<0,2,3> fe_p0;
//     FE_P0_XFEM<2,3> fe_p0_xfem(&fe_p0,{func});
//     FEValues<2,3> fe_values_p(map, *qxfem, fe_p0_xfem, update_values | update_JxW_values | update_inverse_jacobians | update_quadrature_points);
//     
//     fe_values_p.reinit(ele);
//     
//     cout << "Side values:" << endl;
//     QGauss<1> side_quad(1);
//     FESideValues<2,3> fv_side(map, side_quad, fe_rt_xfem, update_values | update_JxW_values
//                                                           | update_quadrature_points | update_normal_vectors);
//     
//     print_fv_side(ele,fv_side, 0);
//     print_fv_side(ele,fv_side, 1);
//     print_fv_side(ele,fv_side, 2);
// //     fv_side.reinit(ele, 0);
// //     cout << fv_side.n_points() << endl;
// //     cout << fv_side.n_dofs() << endl;
// //     cout << fv_side.shape_vector(0,0) << endl;
// //     cout << fv_side.normal_vector(0) << endl;
// //     cout << side_quad.weight(0) << endl;
// //     cout << ele->side(0)->measure() << endl;
// //     cout << arma::dot(fv_side.shape_vector(0,0),fv_side.normal_vector(0))*ele->side(0)->measure() << endl;
// //     cout << arma::dot(fv_side.shape_vector(1,0),fv_side.normal_vector(0))*ele->side(0)->measure() << endl;
// //     cout << arma::dot(fv_side.shape_vector(2,0),fv_side.normal_vector(0))*ele->side(0)->measure() << endl;
// //     cout << arma::dot(fv_side.shape_vector(3,0),fv_side.normal_vector(0))*ele->side(0)->measure() << endl;
// //     cout << arma::dot(fv_side.shape_vector(4,0),fv_side.normal_vector(0))*ele->side(0)->measure() << endl;
// //     cout << arma::dot(fv_side.shape_vector(5,0),fv_side.normal_vector(0))*ele->side(0)->measure() << endl;
//     
// }



// // simplest mesh
// string ref_element_mesh2 = R"CODE(
// $MeshFormat
// 2.2 0 8
// $EndMeshFormat
// $Nodes
// 3
// 1 2.5 1.875 0
// 2 2.5 2.5 0
// 3 1.875 2.5 0
// $EndNodes
// $Elements
// 1
// 1 2 2 39 40 1 3 2
// $EndElements
// )CODE";

// 1 0 0 0
// 2 2 0 0
// 3 0 2 0

// void enriched_side_edge(ElementFullIter ele, FE_RT0_XFEM<2,3>* fe_rt_xfem_, unsigned int local_side){
//         double val;
//         DBGVAR(local_side);
//         MappingP1<2,3> map;
//         
//         QGauss<1> auxq(1);
//         auto fv_side = std::make_shared<FESideValues<2,3>>(map, auxq, *fe_rt_xfem_, update_normal_vectors);
//         fv_side->reinit(ele, local_side);
//         
//         const unsigned int qsize = 100;
//         QXFEM<2,3> qside_xfem(QMidpoint(qsize), local_side, *ele->permutation_idx_); // mapped side quadrature to 2d coords
//         DBGVAR(*ele->permutation_idx_);
//         for(unsigned int q=0; q < qside_xfem.size(); q++){   // map to real coords
//             arma::vec real_point = map.project_unit_to_real(RefElement<2>::local_to_bary(qside_xfem.point(q)),map.element_map(*ele));
//             qside_xfem.set_real_point(q,real_point);
//         }
//         
//         auto fv_xfem = std::make_shared<FEValues<2,3>>(map, qside_xfem, *fe_rt_xfem_, update_values);
//         fv_xfem->reinit(ele);
//         
//         for(unsigned int j=fe_rt_xfem_->n_regular_dofs(); j<fe_rt_xfem_->n_dofs(); j++){
//             double sum_val = 0;
//             double side_measure = ele->side(local_side)->measure();
//             for(unsigned int q=0; q < qside_xfem.size(); q++){
// //                 auto qp = qside_xfem.real_point(q);
// //                 cout << qp(0) << " " << qp(1) << " " << qp(2) << "\n";
// //                 auto fv = fv_xfem->shape_vector(j,q);
// //                 cout << fv(0) << " " << fv(1) << " " << fv(2) << "\n";
//                 val = arma::dot(fv_xfem->shape_vector(j,q),fv_side->normal_vector(0))
//                       // this makes JxW on the triangle side:
//                       * qside_xfem.weight(q)
//                       * side_measure;
//                       
//                 sum_val += val;
//             }
//             DBGVAR(sum_val);
//         }           
//     }
//     
// TEST(fe_xfem, fe_rt_xfem_edges) {
// 
//     QXFEMFactory qfactory;
//     
//     // read mesh - simplset cube from test1
//     Mesh* mesh = mesh_constructor();
//     stringstream in(ref_element_mesh2.c_str());
//     mesh->read_gmsh_from_stream(in);
//     ElementFullIter ele = mesh->element(0);
//     Point n = arma::cross(ele->node[1]->point() - ele->node[0]->point(),
//                           ele->node[2]->point() - ele->node[0]->point());
//     
// //     auto func = std::make_shared<Singularity0D<3>>(arma::vec({0.3,0.3,2}),0.1,arma::vec({0,0,1}),n);
//     auto func = std::make_shared<Singularity0D>(arma::vec({3.3,3.3,0}),0.02,arma::vec({0,0,1}),n, 100);
//     shared_ptr<QXFEM<2,3>> qxfem = qfactory.create_singular({func},ele);
//     
//     MappingP1<2,3> map;
//     
//     FE_RT0<2,3> fe_rt;
//     FE_RT0_XFEM<2,3> fe_rt_xfem(&fe_rt,{func});
//     
//     for(unsigned int j=0; j<3; j++){
//         enriched_side_edge(ele, &fe_rt_xfem, j);
//     }
// }
