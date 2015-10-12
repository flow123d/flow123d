/*
 * intersection_area_test.cpp
 *
 *  Created on: Nov 14, 2014
 *      Author: VF
 */
#define TEST_USE_MPI
#include <flow_gtest_mpi.hh>

#include "system/system.hh"
#include "system/sys_profiler.hh"
#include "system/file_path.hh"
#include "mesh/mesh.h"

#include "mesh/ngh/include/point.h"
#include "mesh/ngh/include/intersection.h"

#include "intersection/inspectelements.h"

using namespace std;
using namespace computeintersection;


TEST(area_intersections, all) {
    
	xprintf(Msg, "==============\n");
	xprintf(Msg, "Computing polygon area by NEW algorithm\n");
	xprintf(Msg, "==============\n");
	double area1, area2;

    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
	string path = "/mesh/site/triangle_tetrahedron.msh";
	FilePath mesh_file(path, FilePath::input_file);
	Profiler::initialize();
	Mesh mesh;
    ifstream ifs(string(mesh_file).c_str());
    mesh.read_gmsh_from_stream(ifs);


	InspectElements ie(&mesh);
// 	ie.print(0);
// 	ie.print(1);
    ie.compute_intersections<2,3>();

	area1 = ie.polygonArea();
	xprintf(Msg, "Polygon area successfully computed\n\n");

	xprintf(Msg, "==============\n");
	xprintf(Msg, "Computing polygon area by NGH algorithm\n");
	xprintf(Msg, "==============\n");

	TTriangle ttr;
	TTetrahedron tte;
	TIntersectionType it;

	FOR_ELEMENTS(&mesh, elm) {
		if (elm->dim() == 2) {
		ttr.SetPoints(TPoint(elm->node[0]->point()(0), elm->node[0]->point()(1), elm->node[0]->point()(2)),
					 TPoint(elm->node[1]->point()(0), elm->node[1]->point()(1), elm->node[1]->point()(2)),
					 TPoint(elm->node[2]->point()(0), elm->node[2]->point()(1), elm->node[2]->point()(2)) );
		}else if(elm->dim() == 3){
		tte.SetPoints(TPoint(elm->node[0]->point()(0), elm->node[0]->point()(1), elm->node[0]->point()(2)),
					 TPoint(elm->node[1]->point()(0), elm->node[1]->point()(1), elm->node[1]->point()(2)),
					 TPoint(elm->node[2]->point()(0), elm->node[2]->point()(1), elm->node[2]->point()(2)),
					 TPoint(elm->node[3]->point()(0), elm->node[3]->point()(1), elm->node[3]->point()(2)));
		}
	}

	GetIntersection(ttr, tte, it, area2);
	xprintf(Msg, "Obsah novy algoritmus: %f, obsah ngh: %f\n", area1, area2);
    
    Profiler::uninitialize();
	xprintf(Msg, "Test complete!\n");
}




