#define TEST_USE_MPI
#include <flow_gtest_mpi.hh>
//#define Flow123d_DEBUG
#include "system/system.hh"
#include "system/sys_profiler.hh"
#include "system/file_path.hh"
#include <array>
#include "mesh/msh_gmshreader.h"
#include "intersection/inspectelements.h"

using namespace std;
using namespace computeintersection;


TEST(intersections, all) {

	cout << "===============" << endl;
	FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
	FilePath mesh_file("mesh/site/TEST.msh", FilePath::input_file);

	Profiler::initialize();

	Mesh mesh;
	ifstream ifs(string(mesh_file).c_str());
	mesh.read_gmsh_from_stream(ifs);

	cout << "Síť načtena!" << endl;
	cout << "Probíhá výpočet průniku" << endl;

	InspectElements ie(&mesh);
	{ START_TIMER("Vypocet pruniku");

	//ie.ComputeIntersections23();
	ie.compute_intersections<2,3>();
	END_TIMER("Vypocet pruniku");}
	//ie.print(0);
	//ie.print(1);
	ie.print_mesh_to_file("pokus");

	double obsah = ie.polygonArea();
	xprintf(Msg,"Obsah polygonu: %f\n", obsah);

	Profiler::instance()->output(MPI_COMM_WORLD, cout);
	Profiler::uninitialize();

	//count << pa.c_str() << endl;
	xprintf(Msg, "Test complete! ");
}




