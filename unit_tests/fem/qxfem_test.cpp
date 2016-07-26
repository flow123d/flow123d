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
1 2 2 39 40 1 3 2
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
    // add this to splot: 'q_points.dat' using 1:2:3 with points lc rgb 'green' title 'q_points'   
}



class TestQXFEMFactory : public testing::Test, public QXFEMFactory<2,3> {
public:
    TestQXFEMFactory(unsigned int max_level = 10) : QXFEMFactory<2,3>(max_level) {}

    void test_distance(const Singularity0D<3>& sing, AuxSimplex& s, const Point& u)
    {
        double computed_distance_sqr = -1.0;
        double max_h = 0;
        simplex_sigularity_intersection(sing,s,computed_distance_sqr, max_h);
        
//         DBGMSG("maxh %f crit %f\n",max_h, square_refinement_criteria_factor_ * computed_distance_sqr);
        
        double accurate_distance_sqr = arma::dot(u,u);
//         DBGMSG("%f %f %e\n",computed_distance_sqr, accurate_distance_sqr, computed_distance_sqr-accurate_distance_sqr);
        EXPECT_NEAR(accurate_distance_sqr, computed_distance_sqr, 1e-14);
        
        
//         string dir_name = string(UNIT_TESTS_SRC_DIR) + "/fem/qxfem_output/";
//         std::ofstream q_points_file;
//         q_points_file.open (dir_name + "distance.dat");
//         if (q_points_file.is_open())
//         {
//             for(const Space<3>::Point &p : s.nodes)
//                 q_points_file << p[0] << " " << p[1] << " " << p[2] << "\n";
//             
//             q_points_file << s.nodes[0][0] << " " << s.nodes[0][1] << " " << s.nodes[0][2] << "\n";
//             
//             q_points_file << sing.center()[0] << " " << sing.center()[1] << " " << sing.center()[2] << "\n";
//         }
//         else
//         {
//             std::cout << "Coud not write refinement for gnuplot.\n";
//         }
//         q_points_file.close();
    }
};


TEST_F(TestQXFEMFactory, distance) {
    // prepare simplex
    AuxSimplex s;
    std::vector<Point> nodes = {{2,3,0}, {8,1,3}, {7,8,6}};
    s.nodes.resize(nodes.size());
    
    const std::vector<std::vector<unsigned int>> permutations_triangle = {
    {0,1,2},
    {1,0,2},
    {1,2,0},
    {0,2,1},
    {2,0,1},
    {2,1,0}};
    
    for(unsigned int p=0; p<permutations_triangle.size(); p++) {        
        for(unsigned int i=0; i<nodes.size(); i++)
        {
            s.nodes[i] = nodes[permutations_triangle[p][i]];
        }
    
        Point v0 = s.nodes[1] - s.nodes[0],   // 0. edge of triangle
            v1 = s.nodes[2] - s.nodes[1],   // 1. edge of triangle
            v2 = s.nodes[0] - s.nodes[2];   // 2. edge of triangle
        
        // normals
        Point nn = arma::cross(v0,v1);
        Point n0 = arma::cross(v0, nn); n0 = n0/arma::norm(n0,2);
        Point n1 = arma::cross(v1, nn); n1 = n1/arma::norm(n1,2);
        Point n2 = arma::cross(v2, nn); n2 = n2/arma::norm(n2,2);

        // length factor
        double f = 1;
        // singularity radius
        double radius = 0.01;
        
        Point m,u,c;
        
        // Singularity v0
        m = s.nodes[0] + 0.5*v0;
        u = f*n0;
        c = m + u;

    //     DBGMSG("k0 = %f\n",arma::dot(m-s.nodes[0], m-s.nodes[0]));
        test_distance(Singularity0D<3>(c,radius), s, u);
        
        // Singularity v1
        m = s.nodes[1] + 0.33*v1;
        u = f*n1;
        c = m + u;

    //     DBGMSG("k2 = %f\n",arma::dot(m-s.nodes[2], m-s.nodes[2]));
        test_distance(Singularity0D<3>(c,radius), s, u);
        
        // Singularity v2
        m = s.nodes[2] + 0.6*v2;
        u = f*n2;
        c = m + u;

    //     DBGMSG("k1 = %f\n",arma::dot(m-s.nodes[1], m-s.nodes[1]));
        test_distance(Singularity0D<3>(c,radius), s, u);
        
        // Singularity node0
        m = s.nodes[0];
        u = -f*n1;
        c = m + u;

        test_distance(Singularity0D<3>(c,radius), s, u);
        
        // Singularity node1
        m = s.nodes[1];
        u = -f*n2;
        c = m + u;

        test_distance(Singularity0D<3>(c,radius), s, u);
        
        // Singularity node2
        m = s.nodes[2];
        u = -f*n0;
        c = m + u;

        test_distance(Singularity0D<3>(c,radius), s, u);
    }
}