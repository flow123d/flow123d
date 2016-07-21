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

#include "fem/singularity.hh"

#include "quadrature/qxfem.hh"
#include "quadrature/qxfem_factory.hh"

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
1 2 2 39 40 1 2 3
$EndElements
)CODE";


TEST(qxfem, singularity) {
    
    // read mesh - simplset cube from test1
    Mesh* mesh = new Mesh();
    stringstream in(ref_element_mesh.c_str());
    mesh->read_gmsh_from_stream(in);
    ElementFullIter ele = mesh->element(0);
 
    Space<3>::Point center({1,2,2});
    double radius = 0.1;
    unsigned int n_qpoints = 100;
    
    Singularity0D<3> func(center,radius);
    func.evaluate_q_points(n_qpoints, ele);
    
    EXPECT_EQ(n_qpoints, func.q_points().size());
    for(const Space<3>::Point &p : func.q_points()){
        double dist = arma::norm(p-center, 2);
        EXPECT_NEAR(radius, dist, 1e-15);
    }
    
//     string dir_name = string(UNIT_TESTS_SRC_DIR) + "/fem/qxfem_output/";
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
    
    Space<2>::Point center2d({1,2});
    Singularity0D<2> func2d(center2d,radius);
    func2d.evaluate_q_points(n_qpoints);
    EXPECT_EQ(n_qpoints, func2d.q_points().size());
    for(const Space<2>::Point &p : func2d.q_points()){
        double dist = arma::norm(p-center2d, 2);
        EXPECT_NEAR(radius, dist, 1e-15);
    }
}


TEST(qxfem, qxfem_factory) {

    QXFEMFactory<2,3> qfactory;
    
    // read mesh - simplset cube from test1
    Mesh* mesh = new Mesh();
    stringstream in(ref_element_mesh.c_str());
    mesh->read_gmsh_from_stream(in);
    ElementFullIter ele = mesh->element(0);
 
    
    Singularity0D<3> func({1,2,2},0.1);
    shared_ptr<QXFEM<2,3>> qxfem2 = qfactory.create_singular({func},ele);
    
//     string dir_name = string(UNIT_TESTS_SRC_DIR) + "/fem/qxfem_output/";
//     qfactory.gnuplot_refinement(ele, dir_name, *qxfem2, {func});
    
    
    double sum=0;
    for(unsigned int q=0; q<qxfem2-> size(); q++) sum += qxfem2->weight(q);
    EXPECT_NEAR(sum,ele->measure()-M_PI*0.01,3e-5);
//     cout << setprecision(15) << "sum weigths: " << sum - (ele->measure()-M_PI*0.01) << endl;
    
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
    // add this to splot: 'q_points.dat' using 1:2:3 with points lc rgb 'green' title 'q_points' ,\
    
}