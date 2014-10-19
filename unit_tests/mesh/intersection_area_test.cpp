/*
 * intersection_area_test.cpp
 *
 *  Created on: Nov 14, 2014
 *      Author: VF
 */
#include <flow_gtest.hh>
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


TEST(area_intersections, all) {

	xprintf(Msg, "==============\n");
	xprintf(Msg, "Computing polygon area by NEW algorithm\n");
	xprintf(Msg, "==============\n");
	double area1, area2;

	string path = (string)UNIT_TESTS_SRC_DIR + "/mesh/intersection_area_test.msh";
	FilePath mesh_file(path, FilePath::input_file);
	Profiler::initialize();
	Mesh mesh;
	GmshMeshReader reader(mesh_file);
	reader.read_mesh(&mesh);

	InspectElements ie(&mesh);
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
	xprintf(Msg, "Test complete!\n");
}




