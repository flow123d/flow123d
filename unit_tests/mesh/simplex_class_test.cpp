#define TEST_USE_MPI
#include <flow_gtest_mpi.hh>
#include "system/system.hh"
//#include "system/sys_profiler.hh"
//#include "system/file_path.hh"
#include <array>
#include "mesh/msh_gmshreader.h"

#include "mesh/ngh/include/point.h"
#include "mesh/ngh/include/intersection.h"

#include "intersection/inspectelements.h"

using namespace std;
using namespace computeintersection;

TEST(simplex, all) {

	//FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");



	arma::vec3 Point0;Point0[0] = 1;Point0[1] = 2;Point0[2] = 3;
	arma::vec3 Point1;Point1[0] = 2;Point1[1] = 3;Point1[2] = 4;
	arma::vec3 Point2;Point2[0] = 3;Point2[1] = 4;Point2[2] = 1;
	arma::vec3 Point3;Point3[0] = 4;Point3[1] = 1;Point3[2] = 2;

	arma::vec3 PointA;PointA[0] = 5;PointA[1] = 6;PointA[2] = 7;
	arma::vec3 PointB;PointB[0] = 6;PointB[1] = 7;PointB[2] = 5;
	arma::vec3 PointC;PointC[0] = 7;PointC[1] = 5;PointC[2] = 6;


	arma::vec3 *pole_pp[] = {&Point0,&Point1,&Point2,&Point3};
	arma::vec3 *pole_p[] = {&PointA, &PointB, &PointC};

	Simplex<3> ss(pole_pp);
	Simplex<2> s(pole_p);


	FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

	Profiler::initialize();

	FilePath mesh_file("mesh/site/triangle_tetrahedron12.msh", FilePath::input_file);

	Mesh mesh;
	ifstream ifs(string(mesh_file).c_str());
	mesh.read_gmsh_from_stream(ifs);


	InspectElements ie(&mesh);
	//ie.ComputeIntersections23();
	//ie.print_mesh_to_file("neco");
	/*ComputeIntersection<Simplex<2>, Simplex<3>> pp(s,ss);
	pp.init();

	pp.toStringPluckerCoordinatesTree();
*/
	Profiler::uninitialize();

	ie.compute_intersections<2,3>();
	ie.compute_intersections<1,2>();
	ie.compute_intersections<0,94>();
	ie.compute_intersections<8,5>();

	cout << "eps: " << 64*numeric_limits<double>::epsilon()<< endl;

	xprintf(Msg, "Test complete!\n");
}




