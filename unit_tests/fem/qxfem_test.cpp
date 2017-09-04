/*
 * 
 * 
 */

#define TEST_USE_PETSC
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>

#include "arma_expect.hh"

#include "mesh/mesh.h"
#include "mesh_constructor.hh"

#include "fem/singularity.hh"

#include "quadrature/qxfem.hh"
#include "quadrature/qxfem_factory.hh"




typedef Space<3>::Point Point;

TEST(qxfem, qxfem_factory) {

    QXFEMFactory qfactory(12);
    
    // setup FilePath directories
    FilePath::set_io_dirs(".",FilePath::get_absolute_working_dir(),"",".");
    string dir_name = string(UNIT_TESTS_SRC_DIR) + "/fem/qxfem_output/";
    // read mesh
    Mesh *mesh = mesh_full_constructor("{mesh_file=\"triangle.msh\"}");
    
    ElementFullIter ele = mesh->element(0);
    Point n = arma::cross(ele->node[1]->point() - ele->node[0]->point(),
                          ele->node[2]->point() - ele->node[0]->point());
    
    auto func = std::make_shared<Singularity<0>>(arma::vec({1,2,2}),0.1,arma::vec({0,0,1}),n, 100);
    shared_ptr<QXFEM<2,3>> qxfem = qfactory.create_singular({func},ele);
    
//     string dir_name = string(UNIT_TESTS_SRC_DIR) + "/fem/qxfem_output/";
//     qfactory.gnuplot_refinement(ele, dir_name, *qxfem);
    
//     std::ofstream q_points_file;
//     q_points_file.open (dir_name + "unit_q_points.dat");
//     if (q_points_file.is_open()) 
//     {
//         for(const Space<2>::Point &p : qxfem->get_points())
//             q_points_file << p[0] << " " << p[1] << " " << 0 << "\n";
//     }
//     else 
//     { 
//         MessageOut() << "Coud not write refinement for gnuplot.\n";
//     }
//     q_points_file.close();
    
    double sum=0;
    for(unsigned int q=0; q<qxfem-> size(); q++) sum += qxfem->weight(q);
//     MessageOut() << setprecision(15) << "sum: " << sum << "\n";
//     MessageOut() << setprecision(15) << "Tmeasure: " << ele->measure() << "\n";
    
    double exact_sum = (ele->measure()-func->geometry().volume()) / (2*ele->measure());
//     MessageOut() << setprecision(15) << "exact_sum: " << exact_sum << "\n";
    MessageOut() << setprecision(15) << "sum weigths diff: " << sum - exact_sum << "\n";
    EXPECT_NEAR(sum,exact_sum,1e-7);
    
//     std::ofstream q_points_file;
//     q_points_file.open (dir_name + "q_points.dat");
//     if (q_points_file.is_open()) 
//     {
//         for(const Singularity<0>::Point &p : func.q_points())
//         q_points_file << p[0] << " " << p[1] << " " << p[2] << "\n";
//     }
//     else 
//     { 
//         WarningOut() << "Coud not write refinement for gnuplot.\n";
//     }
//     q_points_file.close();
    // add this to splot: 'q_points.dat' using 1:2:3 with points lc rgb 'green' title 'q_points'   
}



class TestQXFEMFactory : public testing::Test, public QXFEMFactory {
public:
    TestQXFEMFactory(unsigned int max_level = 10) : QXFEMFactory(max_level) {}

    void test_distance(const Singularity<0>& sing, AuxSimplex& s, const Point& u)
    {
        double computed_distance_sqr = -1.0;
//         double max_h = s.compute_max_h<2>();
        sigularity0D_distance(sing,s,computed_distance_sqr);
        
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
//             WarningOut() << "Coud not write refinement for gnuplot.\n";
//         }
//         q_points_file.close();
    }
};


TEST_F(TestQXFEMFactory, distance) {
    // prepare simplex
    AuxSimplex s;
    std::vector<Point> nodes = {{2,3,0}, {8,1,3}, {7,8,6}};
    s.nodes.resize(nodes.size());
    unsigned int n_qpoints = 100;
    
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
        Space<3>::Point dv({2,2,2});
        
        Point m,u,c;
        
        // Singularity v0
        m = s.nodes[0] + 0.5*v0;
        u = f*n0;
        c = m + u;

    //     DBGMSG("k0 = %f\n",arma::dot(m-s.nodes[0], m-s.nodes[0]));
        test_distance(Singularity<0>(c,radius,dv,nn,n_qpoints), s, u);
        
        // Singularity v1
        m = s.nodes[1] + 0.33*v1;
        u = f*n1;
        c = m + u;

    //     DBGMSG("k2 = %f\n",arma::dot(m-s.nodes[2], m-s.nodes[2]));
        test_distance(Singularity<0>(c,radius,dv,nn,n_qpoints), s, u);
        
        // Singularity v2
        m = s.nodes[2] + 0.6*v2;
        u = f*n2;
        c = m + u;

    //     DBGMSG("k1 = %f\n",arma::dot(m-s.nodes[1], m-s.nodes[1]));
        test_distance(Singularity<0>(c,radius,dv,nn,n_qpoints), s, u);
        
        // Singularity node0
        m = s.nodes[0];
        u = -f*n1;
        c = m + u;

        test_distance(Singularity<0>(c,radius,dv,nn,n_qpoints), s, u);
        
        // Singularity node1
        m = s.nodes[1];
        u = -f*n2;
        c = m + u;

        test_distance(Singularity<0>(c,radius,dv,nn,n_qpoints), s, u);
        
        // Singularity node2
        m = s.nodes[2];
        u = -f*n0;
        c = m + u;

        test_distance(Singularity<0>(c,radius,dv,nn,n_qpoints), s, u);
    }
}



TEST(qxfem, qxfem_factory_two) {

    QXFEMFactory qfactory(12);
    
    // setup FilePath directories
    FilePath::set_io_dirs(".",FilePath::get_absolute_working_dir(),"",".");
    string dir_name = string(UNIT_TESTS_SRC_DIR) + "/fem/qxfem_output/";
    // read mesh
    Mesh *mesh = mesh_full_constructor("{mesh_file=\"triangle.msh\"}");
    
    ElementFullIter ele = mesh->element(0);
    Point n = arma::cross(ele->node[1]->point() - ele->node[0]->point(),
                          ele->node[2]->point() - ele->node[0]->point());
    unsigned int n_qpoints = 100;
    
    auto func1 = std::make_shared<Singularity<0>>(arma::vec({1,2,2}),0.1,arma::vec({0,0,1}),n,n_qpoints);
    auto func2 = std::make_shared<Singularity<0>>(arma::vec({2,1,2}),0.05,arma::vec({0,0,1}),n,n_qpoints);
    shared_ptr<QXFEM<2,3>> qxfem = qfactory.create_singular({func1, func2},ele);
    
//     string dir_name = string(UNIT_TESTS_SRC_DIR) + "/fem/qxfem_output_two/";
//     qfactory.gnuplot_refinement(ele, dir_name, *qxfem, {func1, func2});
//     
//     std::ofstream q_points_file;
//     q_points_file.open (dir_name + "unit_q_points.dat");
//     if (q_points_file.is_open()) 
//     {
//         for(const Space<2>::Point &p : qxfem->get_points())
//             q_points_file << p[0] << " " << p[1] << " " << 0 << "\n";
//     }
//     else 
//     { 
//         MessageOut() << "Coud not write refinement for gnuplot.\n";
//     }
//     q_points_file.close();
    
    double sum=0;
    for(unsigned int q=0; q<qxfem-> size(); q++) sum += qxfem->weight(q);
//     MessageOut() << setprecision(15) << "sum: " << sum << "\n";
//     MessageOut() << setprecision(15) << "Tmeasure: " << ele->measure() << "\n";
    
    double exact_sum = (ele->measure()-func1->geometry().volume()-func2->geometry().volume()) / (2*ele->measure());
//     MessageOut() << setprecision(15) << "exact_sum: " << exact_sum << "\n";
    MessageOut() << setprecision(15) << "sum weigths diff: " << sum - exact_sum << "\n";
    EXPECT_NEAR(sum,exact_sum,1e-7);
    
    
//     func.evaluate_q_points(100, ele);
//     std::ofstream q_points_file;
//     q_points_file.open (dir_name + "q_points.dat");
//     if (q_points_file.is_open()) 
//     {
//         for(const Singularity<0>::Point &p : func.q_points())
//         q_points_file << p[0] << " " << p[1] << " " << p[2] << "\n";
//     }
//     else 
//     { 
//         WarningOut() << "Coud not write refinement for gnuplot.\n";
//     }
//     q_points_file.close();
    // add this to splot: 'q_points.dat' using 1:2:3 with points lc rgb 'green' title 'q_points'   
}

















TEST(qxfem, qxfem_factory_3d) {

    // setup FilePath directories
    FilePath::set_io_dirs(".",FilePath::get_absolute_working_dir(),"",".");
    string dir_name = string(UNIT_TESTS_SRC_DIR) + "/fem/qxfem_output/";
    // read mesh
    Mesh *mesh = mesh_full_constructor("{mesh_file=\"tetrahedron.msh\"}");
    
    ElementFullIter ele = mesh->element(0);
    
    unsigned int n = 50, m = 20;
    
//     auto func = std::make_shared<Singularity<1>>(arma::vec({0.75,0.5,0}), arma::vec({0.75,0.5,3}),0.2);
    auto func = std::make_shared<Singularity<1>>(arma::vec({0.5,0.5,0}), arma::vec({0.5,0.5,2}),0.05,n,m);
    
//     for(unsigned int i=4; i<10; i++)
    unsigned int i = 9;
    {
        QXFEMFactory qfactory(i);
        shared_ptr<QXFEM<3,3>> qxfem = qfactory.create_singular({func},ele);
//         qfactory.gnuplot_refinement<3>(ele, dir_name, *qxfem);
        
    //     std::ofstream q_points_file;
    //     q_points_file.open (dir_name + "unit_q_points.dat");
    //     if (q_points_file.is_open()) 
    //     {
    //         for(const Space<2>::Point &p : qxfem->get_points())
    //             q_points_file << p[0] << " " << p[1] << " " << 0 << "\n";
    //     }
    //     else 
    //     { 
    //         MessageOut() << "Coud not write refinement for gnuplot.\n";
    //     }
    //     q_points_file.close();
        
        double sum=0;
        for(unsigned int q=0; q<qxfem->size(); q++) sum += qxfem->weight(q);
//         MessageOut() << setprecision(15) << "sum: " << sum << "\n";
//         MessageOut() << setprecision(15) << "Tmeasure: " << ele->measure() << "\n";
//         MessageOut() << setprecision(15) << "Cylinder: " << func->geometry().volume() << "\n";
        
        double exact_sum = (ele->measure()-func->geometry().volume()*0.5) / (6*ele->measure());
        MessageOut() << setprecision(15) << "exact_sum: " << exact_sum << "\n";
        MessageOut() << setprecision(15) << "sum weigths diff: " << sum - exact_sum << "\n";
        EXPECT_NEAR(sum,exact_sum,1e-6);
    }
    
//     std::ofstream q_points_file;
//     q_points_file.open (dir_name + "q_points.dat");
//     if (q_points_file.is_open()) 
//     {
//         for(const Singularity<1>::Point &p : func->q_points())
//         q_points_file << p[0] << " " << p[1] << " " << p[2] << "\n";
//     }
//     else 
//     { 
//         WarningOut() << "Coud not write refinement for gnuplot.\n";
//     }
//     q_points_file.close();
    // add this to splot: 'q_points.dat' using 1:2:3 with points lc rgb 'green' title 'q_points'
}

TEST(qxfem, qxfem_factory_3d_side) {

    // setup FilePath directories
    FilePath::set_io_dirs(".",FilePath::get_absolute_working_dir(),"",".");
    string dir_name = string(UNIT_TESTS_SRC_DIR) + "/fem/qxfem_output/";
    // read mesh
    Mesh *mesh = mesh_full_constructor("{mesh_file=\"tetrahedron.msh\"}");
    
    ElementFullIter ele = mesh->element(0);
    
    unsigned int n = 50, m = 20;
    
//     auto func = std::make_shared<Singularity<1>>(arma::vec({0.75,0.5,0}), arma::vec({0.75,0.5,3}),0.2);
    auto func = std::make_shared<Singularity<1>>(arma::vec({0.5,0.2,0}), arma::vec({0.2,0.5,2}),0.05,n,m);
    
//     for(unsigned int i=6; i<7; i++)
    unsigned int sid = 3;
    unsigned int i = 12;
    {
        QXFEMFactory qfactory(i);
        shared_ptr<QXFEM<3,3>> qxfem = qfactory.create_side_singular({func},ele, sid);
        
        qfactory.gnuplot_refinement<3>(ele, dir_name, *qxfem);
        
        //comparison with 2d - we need only the are of ellipse
        auto nodes = RefElement<3>::interact(Interaction<0,2>(sid));
        Point n = arma::cross(ele->node[nodes[1]]->point() - ele->node[nodes[0]]->point(),
                              ele->node[nodes[2]]->point() - ele->node[nodes[0]]->point());
        const CylinderGeometry& geom = func->geometry_cylinder();
        
        CircleEllipseProjection geom_ellipse(arma::vec({0,0,0}),geom.radius(),geom.direction_vector(),n);

        
//         std::ofstream q_points_file;
//         q_points_file.open (dir_name + "unit_q_points.dat");
//         if (q_points_file.is_open()) 
//         {
//             for(const Space<2>::Point &p : qxfem->get_points())
//                 q_points_file << p[0] << " " << p[1] << " " << 0 << "\n";
//         }
//         else 
//         { 
//             MessageOut() << "Coud not write refinement for gnuplot.\n";
//         }
//         q_points_file.close();
        
        double sum=0;
        for(unsigned int q=0; q<qxfem->size(); q++) sum += qxfem->weight(q);
        MessageOut() << setprecision(15) << "sum: " << sum << "\n";
//         MessageOut() << setprecision(15) << "Tmeasure: " << ele->side(sid)->measure() << "\n";
//         MessageOut() << setprecision(15) << "Cylinder: " << func->geometry().volume() << "\n";
        
        double exact_sum = (ele->side(sid)->measure()-geom_ellipse.ellipse_area()) / (2*ele->side(sid)->measure());
        MessageOut() << setprecision(15) << "exact_sum: " << exact_sum << "\n";
        MessageOut() << setprecision(15) << "sum weigths diff: " << sum - exact_sum << "\n";
        EXPECT_NEAR(sum,exact_sum,1e-8);
    }
    
//     func->evaluate_q_points(50, 20);
//     std::ofstream q_points_file;
//     q_points_file.open (dir_name + "q_points.dat");
//     if (q_points_file.is_open()) 
//     {
//         for(const Singularity<1>::Point &p : func->q_points())
//         q_points_file << p[0] << " " << p[1] << " " << p[2] << "\n";
//     }
//     else 
//     { 
//         WarningOut() << "Coud not write refinement for gnuplot.\n";
//     }
//     q_points_file.close();
    // add this to splot: 'q_points.dat' using 1:2:3 with points lc rgb 'green' title 'q_points'
}


TEST(qxfem, qxfem_factory_2d_side) {

    // setup FilePath directories
    FilePath::set_io_dirs(".",FilePath::get_absolute_working_dir(),"",".");
    string dir_name = string(UNIT_TESTS_SRC_DIR) + "/fem/qxfem_output/";
    // read mesh
    Mesh *mesh = mesh_full_constructor("{mesh_file=\"triangle.msh\"}");
    
    ElementFullIter ele = mesh->element(0);
    Point n = arma::cross(ele->node[1]->point() - ele->node[0]->point(),
                          ele->node[2]->point() - ele->node[0]->point());
    
    unsigned int n_qpoints = 100;
    auto func = std::make_shared<Singularity<0>>(arma::vec({1.4,2.9,3}),0.1,arma::vec({1.4,2.9,0}),n,n_qpoints);
    
    unsigned int sid = 2;
    unsigned int i = 12;
    QXFEMFactory qfactory(i);
    shared_ptr<QXFEM<2,3>> qxfem = qfactory.create_side_singular({func},ele,sid);
    
//     qfactory.gnuplot_refinement<2>(ele, dir_name, *qxfem);
    
//     std::ofstream q_points_file;
//     q_points_file.open (dir_name + "unit_q_points.dat");
//     if (q_points_file.is_open()) 
//     {
//         for(const Space<2>::Point &p : qxfem->get_points())
//             q_points_file << p[0] << " " << p[1] << " " << 0 << "\n";
//     }
//     else 
//     { 
//         MessageOut() << "Coud not write refinement for gnuplot.\n";
//     }
//     q_points_file.close();
    
    double sum=0;
    for(unsigned int q=0; q<qxfem-> size(); q++) sum += qxfem->weight(q);
    MessageOut() << setprecision(15) << "sum: " << sum << "\n";
//     MessageOut() << setprecision(15) << "Tmeasure: " << ele->side(sid)->measure() << "\n";
    
//     double exact_sum = (ele->measure()-func->geometry().ellipse_area()) / (2*ele->measure());
//     MessageOut() << setprecision(15) << "exact_sum: " << exact_sum << "\n";
//     MessageOut() << setprecision(15) << "sum weigths diff: " << sum - exact_sum << "\n";
//     EXPECT_NEAR(sum,exact_sum,1e-7);
    
    
//     func->evaluate_q_points(100);
//     std::ofstream q_points_file;
//     q_points_file.open (dir_name + "q_points.dat");
//     if (q_points_file.is_open()) 
//     {
//         for(const Singularity<0>::Point &p : func->q_points())
//         q_points_file << p[0] << " " << p[1] << " " << p[2] << "\n";
//     }
//     else 
//     { 
//         WarningOut() << "Coud not write refinement for gnuplot.\n";
//     }
//     q_points_file.close();
    // add this to splot: 'q_points.dat' using 1:2:3 with points lc rgb 'green' title 'q_points'   
}