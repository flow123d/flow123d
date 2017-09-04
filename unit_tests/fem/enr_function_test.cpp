/*
 * 
 * 
 */

#define TEST_USE_PETSC
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest.hh>

#include "arma_expect.hh"

#include "mesh/point.hh"
#include "mesh/elements.h"
#include "mesh/mesh.h"
#include "mesh_constructor.hh"

#include "fem/singularity.hh"

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
//         print_point(MessageOut() << setprecision(15),points[1]);
//     }
//     else 
//     { 
//         WarningOut() << "Coud not write for gnuplot.\n";
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



TEST(qxfem, singularity0D) {
    
    //element nodes
    std::vector<Point> node({{0,0,0},
                             {0,6,4},
                             {3,0,2}});
    
    Point n = arma::cross(node[1] - node[0],
                          node[2] - node[0]);
    
    Point center({1,2,2});
    double radius = 0.1;
    Point direction_vector({0,0,1});
    unsigned int n_qpoints = 10;
    
    Singularity<0> func(center, radius, direction_vector, n, n_qpoints);
    const CircleEllipseProjection & geom = func.geometry_ellipse();
    
    EXPECT_EQ(n_qpoints, geom.q_points().size());
    
    std::vector<Point> points = geom.q_points();
    geom.project_to_circle_plane(points);
    
    for(const Point &p : points){
        double dist = geom.distance(p);
        EXPECT_NEAR(radius, dist, 1e-15);
    }
    
    geom.project_to_ellipse_plane(points);
    
    for(unsigned int i=0; i< points.size(); i++){
        EXPECT_ARMA_EQ(geom.q_points()[i], points[i]);
    }
    
//     string dir_name = string(UNIT_TESTS_SRC_DIR) + "/fem/qxfem_output/";
//     std::ofstream q_points_file;
//     q_points_file.open (dir_name + "q_points.dat");
//     if (q_points_file.is_open()) 
//     {
//         for(const Singularity0D::Point &p : geom.q_points())
//             print_point(q_points_file,p);
//     }
//     else 
//     { 
//         WarningOut() << "Coud not write q_points for gnuplot.\n";
//     }
//     q_points_file.close();
// //     add this to splot: 'q_points.dat' using 1:2:3 with points lc rgb 'green' title 'q_points'
}


TEST(qxfem, singularity1D) {
    
    Point a({1,2,2}), b({3,4,5});
    double radius = 0.1;
    unsigned int n = 100, m = 10;
    
    Singularity<1> func(a, b, radius, n, m);
    const CylinderGeometry & geom = func.geometry_cylinder();
    
    EXPECT_EQ(n*m, geom.q_points().size());
    
    for(unsigned int i=0; i< geom.q_points().size(); i++){
        EXPECT_NEAR(func.geometry().distance(geom.q_points()[i]),radius,1e-14);
    }
    
//     string dir_name = string(UNIT_TESTS_SRC_DIR) + "/fem/qxfem_output/";
//     std::ofstream q_points_file;
//     q_points_file.open (dir_name + "s1_q_points.dat");
//     if (q_points_file.is_open()) 
//     {
//         for(const Point &p : func.q_points())
//             print_point(q_points_file,p);
//     }
//     else 
//     { 
//         WarningOut() << "Coud not write refinement for gnuplot.\n";
//     }
//     q_points_file.close();
// //     add this to splot: 's1_q_points.dat' using 1:2:3 with points lc rgb 'green' title 'q_points'
}