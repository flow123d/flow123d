/*
 * 
 * 
 */

#define TEST_USE_PETSC
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>

#include "arma_expect.hh"

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

typedef Space<3>::Point Point;

void print_point(ostream& stream,Point p) {
    stream << p[0] << " " << p[1] << " " << p[2] << "\n";
}

class TestProjection : public testing::Test, public CircleEllipseProjection {
public:
    TestProjection() : CircleEllipseProjection({4,5,6}, 2, {1,2,3}, {5,-2,2}) {}

protected:
    Point c = {4,5,6};
    double r = 2;
    Point u = {1,2,3};
    Point n = {5,-2,2};
};

TEST_F(TestProjection, project_circle_ellipse) {

    EXPECT_EQ(M_PI*r*r, this->circle_area());
    
    u = u / arma::norm(u,2);
    n = n / arma::norm(n,2);
    double cos_a = arma::dot(n,u);
    
    EXPECT_EQ(cos_a, this->cos_a);
    
    EXPECT_EQ(r, this->a_);
    EXPECT_EQ(r/cos_a, this->b_);
    
    Point v = arma::cross(n,u);
    Point w = arma::cross(n,v);
    v = v / arma::norm(v,2);
    w = w / arma::norm(w,2);
    EXPECT_ARMA_EQ(v, this->ea_);
    EXPECT_ARMA_EQ(w, this->eb_);
    
    EXPECT_EQ(M_PI * r * r / cos_a, this->ellipse_area());
    
    std::vector<Point> ell_points, c_points;
    ell_points.push_back(c + ea_ + eb_);
    ell_points.push_back(c + b_*eb_);
    ell_points.push_back(c + a_*ea_);
    ell_points.push_back(c + 0.4*ea_ + 0.3*eb_);
    ell_points.push_back(c + a_*ea_ + 0.0001*eb_);
    
    c_points = ell_points;
    this->project_to_circle_plane(c_points);
    EXPECT_TRUE(point_in_circle(c_points[0]));
    EXPECT_TRUE(point_in_circle(c_points[1]));
    EXPECT_TRUE(point_in_circle(c_points[2]));
    EXPECT_TRUE(point_in_circle(c_points[3]));
    EXPECT_FALSE(point_in_circle(c_points[4]));
    
    this->project_to_ellipse_plane(c_points);
    for(unsigned int i=0; i<ell_points.size(); i++)
        EXPECT_ARMA_EQ(ell_points[i],c_points[i]);

    EXPECT_TRUE(point_in_ellipse(ell_points[0]));
    EXPECT_TRUE(point_in_ellipse(ell_points[1]));
    EXPECT_TRUE(point_in_ellipse(ell_points[2]));
    EXPECT_TRUE(point_in_ellipse(ell_points[3]));
    EXPECT_FALSE(point_in_ellipse(ell_points[4]));
    
    
    
//     Point er = arma::cross(v,u);
//     er = er / arma::norm(er,2);
//     unsigned int count = 100;
//     double phi = 2*M_PI / count;
//     std::vector<Point> circle_points(count), ell_points(count);
//     for(unsigned int i=0; i < count; i++){
//         circle_points[i] = c+ r*sin(i*phi)*ea_+ r*cos(i*phi)*er;
//         ell_points[i] = c + this->a_*cos(i*phi)*this->ea_+ this->b_*sin(i*phi)*this->eb_;
//     }
//     
//     string dir_name = string(UNIT_TESTS_SRC_DIR) + "/fem/qxfem_output/";
//     std::ofstream file;
//     file.open (dir_name + "projection.dat");
//     if (file.is_open())
//     {
//         for(Point& p : circle_points)
//             print_point(file,p);
//         print_point(file,circle_points[0]);
//         
//         file << "\n";
//         for(Point& p : ell_points)
//             print_point(file,p);
//         print_point(file,ell_points[0]);
//         
//         file << "\n";
//         print_point(file,c+2*u);
//         print_point(file,c);
//         
//         file << "\n";
//         print_point(file,c+2*n);
//         print_point(file,c);
//         
//         file << "\n";
//         print_point(file,c + a_*ea_);
//         print_point(file,c);
//         file << "\n";
//         print_point(file,c + b_*eb_);
//         print_point(file,c);
//         
//         file << "\n";
//         print_point(file,c + r*ea_);
//         print_point(file,c);
//         
//         file << "\n";
//         print_point(file,c + r*er);
//         print_point(file,c);
//         
//         
//         file << "\n";
//         print_point(file,points[1]);
//         
//         this->project_to_circle_plane(points);
//         print_point(file,points[1]);
//         print_point(cout << setprecision(15),points[1]);
//     }
//     else 
//     { 
//         std::cout << "Coud not write for gnuplot.\n";
//     }
//     file.close();
}

TEST(CircleEllipseProjection, parallel){
    
    // TEST when direction and normal is parallel
    Point c = {1,2,3};
    double r = 2;
    Point u = {1,2,-1};
    Point n = {2,4,-2};
    CircleEllipseProjection proj(c,r,u,n);
    
    EXPECT_EQ(M_PI * r * r, proj.circle_area());
//     DBGMSG("%e\n",proj.circle_area()- proj.ellipse_area());
    EXPECT_DOUBLE_EQ(proj.circle_area(), proj.ellipse_area());
    
    std::vector<Point> points = {
        {0,0,0},
        {1,2,3},
        {2,2,2},
        {3,2,1}
    };
    auto orig = points;
    
    proj.project_to_circle_plane(points);
    for(unsigned int i=0; i<points.size(); i++)
        EXPECT_ARMA_EQ(orig[i],points[i]);
    
    proj.project_to_ellipse_plane(points);
    for(unsigned int i=0; i<points.size(); i++)
        EXPECT_ARMA_EQ(orig[i],points[i]);
    
    
    // TEST when direction and normal and Z is parallel
    u = {0,0,1};
    n = {0,0,5};
    proj = CircleEllipseProjection(c,r,u,n);
    
    EXPECT_EQ(M_PI * r * r, proj.circle_area());
//     DBGMSG("%e\n",proj.circle_area()- proj.ellipse_area());
    EXPECT_DOUBLE_EQ(proj.circle_area(), proj.ellipse_area());
    
    points = orig;
    
    proj.project_to_circle_plane(points);
    for(unsigned int i=0; i<points.size(); i++)
        EXPECT_ARMA_EQ(orig[i],points[i]);
    
    proj.project_to_ellipse_plane(points);
    for(unsigned int i=0; i<points.size(); i++)
        EXPECT_ARMA_EQ(orig[i],points[i]);
}



TEST(qxfem, singularity) {
    
    // read mesh - simplset cube from test1
    Mesh* mesh = new Mesh();
    stringstream in(ref_element_mesh.c_str());
    mesh->read_gmsh_from_stream(in);
    ElementFullIter ele = mesh->element(0);
    Point n = arma::cross(ele->node[1]->point() - ele->node[0]->point(),
                          ele->node[2]->point() - ele->node[0]->point());
    
    Point center({1,2,2});
    double radius = 0.1;
    Point direction_vector({0,0,1});
    unsigned int n_qpoints = 100;
    
    Singularity0D<3> func(center, radius, direction_vector, n);
    func.evaluate_q_points(n_qpoints);
    
    EXPECT_EQ(n_qpoints, func.q_points().size());
    
    std::vector<Point> points = func.q_points();
    
    func.geometry().project_to_circle_plane(points);
    
    for(const Point &p : points){
        double dist = arma::norm(p-center, 2);
        EXPECT_NEAR(radius, dist, 1e-15);
    }
    
    func.geometry().project_to_ellipse_plane(points);
    for(unsigned int i=0; i< points.size(); i++){
        EXPECT_ARMA_EQ(func.q_points()[i], points[i]);
    }
    
//     string dir_name = string(UNIT_TESTS_SRC_DIR) + "/fem/qxfem_output/";
//     std::ofstream q_points_file;
//     q_points_file.open (dir_name + "q_points.dat");
//     if (q_points_file.is_open()) 
//     {
//         for(const Singularity0D<3>::Point &p : func.q_points())
//             print_point(q_points_file,p);
//     }
//     else 
//     { 
//         std::cout << "Coud not write refinement for gnuplot.\n";
//     }
//     q_points_file.close();
// //     add this to splot: 'q_points.dat' using 1:2:3 with points lc rgb 'green' title 'q_points'
    
    Space<2>::Point center2d({1,2});
    Singularity0D<2> func2d(center2d,radius);
    func2d.evaluate_q_points(n_qpoints);
    EXPECT_EQ(n_qpoints, func2d.q_points().size());
    for(const Space<2>::Point &p : func2d.q_points()){
        double dist = arma::norm(p-center2d, 2);
        EXPECT_NEAR(radius, dist, 1e-15);
    }
}


// void project_circle(const Singularity0D<3>& sing, ElementFullIter ele,
//                     double& a, double& b, double& surface)
// {
//     typedef Space<3>::Point Point;
//     
//     //create first auxsimplex:
//     std::vector<Point> nodes;
//     for(unsigned int i=0; i < ele->n_nodes(); i++)
//         nodes.push_back(ele->node[i]->point());
//     
//     //we suppose counterclockwise node order
//     Point v0 = nodes[1] - nodes[0],   // 0. edge of triangle
//           v1 = nodes[2] - nodes[1];   // 1. edge of triangle
//     if(v0[1]*v1[2] - v0[2]*v1[1] + v0[2]*v1[0] - v0[0]*v1[2] + v0[0]*v1[1] - v0[1]*v1[0] < 0)
//     {
//         DBGMSG("change node order\n");
//         std::swap(nodes[0],nodes[1]);
//     }
//     
//     Point n = arma::cross(v0,v1);   //normal to element
//     Point u = sing.geometry().direction_vector() / arma::norm(sing.geometry().direction_vector(),2);
//     
//     Point v = arma::cross(n,u);
//     Point w = arma::cross(n,v);
//     v = v / arma::norm(v,2) * sing.radius();
//     w = w / arma::norm(w,2);
//     
//     
//     double cos_a = arma::dot(n,u) / arma::norm(n,2);
//     w = sing.radius()/cos_a * w;
//     //double sin_a = std::sqrt(1-cos_a*cos_a);
//     //double tan_a = sin_a / cos_a;
//     
//     (sing.center() + u).print(cout,"direction");
//     (sing.center() + n/arma::norm(n,2)).print(cout,"normal");
//     sing.center().print(cout,"center");
//     (sing.center()+v).print(cout,"a");
//     (sing.center()+w).print(cout,"b");
//     
//     a = arma::norm(v,2);
//     b = arma::norm(w,2);
//     surface = M_PI * a * b;
// }


TEST(qxfem, qxfem_factory) {

    QXFEMFactory<2,3> qfactory(10);
    
    // read mesh - simplset cube from test1
    Mesh* mesh = new Mesh();
    stringstream in(ref_element_mesh.c_str());
    mesh->read_gmsh_from_stream(in);
    ElementFullIter ele = mesh->element(0);
    Point n = arma::cross(ele->node[1]->point() - ele->node[0]->point(),
                          ele->node[2]->point() - ele->node[0]->point());
    
    Singularity0D<3> func({1,2,2},0.1,{0,0,1},n);
    shared_ptr<QXFEM<2,3>> qxfem = qfactory.create_singular({func},ele);
    
    string dir_name = string(UNIT_TESTS_SRC_DIR) + "/fem/qxfem_output/";
    qfactory.gnuplot_refinement(ele, dir_name, *qxfem, {func});
    
    std::ofstream q_points_file;
    q_points_file.open (dir_name + "unit_q_points.dat");
    if (q_points_file.is_open()) 
    {
        for(const Space<2>::Point &p : qxfem->get_points())
            q_points_file << p[0] << " " << p[1] << " " << 0 << "\n";
    }
    else 
    { 
        std::cout << "Coud not write refinement for gnuplot.\n";
    }
    q_points_file.close();
    
    double sum=0;
    for(unsigned int q=0; q<qxfem-> size(); q++) sum += qxfem->weight(q);
    
//     cout << setprecision(15) << "sum weigths diff: " << sum - (ele->measure()-func.geometry().ellipse_area()) << endl;
    EXPECT_NEAR(sum,ele->measure()-func.geometry().ellipse_area(),1e-5);
    
    
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
        Space<3>::Point dv({2,2,2});
        
        Point m,u,c;
        
        // Singularity v0
        m = s.nodes[0] + 0.5*v0;
        u = f*n0;
        c = m + u;

    //     DBGMSG("k0 = %f\n",arma::dot(m-s.nodes[0], m-s.nodes[0]));
        test_distance(Singularity0D<3>(c,radius,dv,nn), s, u);
        
        // Singularity v1
        m = s.nodes[1] + 0.33*v1;
        u = f*n1;
        c = m + u;

    //     DBGMSG("k2 = %f\n",arma::dot(m-s.nodes[2], m-s.nodes[2]));
        test_distance(Singularity0D<3>(c,radius,dv,nn), s, u);
        
        // Singularity v2
        m = s.nodes[2] + 0.6*v2;
        u = f*n2;
        c = m + u;

    //     DBGMSG("k1 = %f\n",arma::dot(m-s.nodes[1], m-s.nodes[1]));
        test_distance(Singularity0D<3>(c,radius,dv,nn), s, u);
        
        // Singularity node0
        m = s.nodes[0];
        u = -f*n1;
        c = m + u;

        test_distance(Singularity0D<3>(c,radius,dv,nn), s, u);
        
        // Singularity node1
        m = s.nodes[1];
        u = -f*n2;
        c = m + u;

        test_distance(Singularity0D<3>(c,radius,dv,nn), s, u);
        
        // Singularity node2
        m = s.nodes[2];
        u = -f*n0;
        c = m + u;

        test_distance(Singularity0D<3>(c,radius,dv,nn), s, u);
    }
}