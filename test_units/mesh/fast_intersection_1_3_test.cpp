/*
 * muj_test.cpp
 *
 *  Created on: 6.2.2013
 *      Author: viktor
 */
#include <gtest/gtest.h>
#include "system/system.hh"
#include "system/file_path.hh"
#include "mesh/mesh.h"
#include "mesh/msh_gmshreader.h"
#include "mesh/bih_tree.hh"
#include "fields/field_interpolated_p0.hh"
#include "system/sys_profiler.hh"

#include "fast_intersection/inspect_elements.h"


using namespace std;
using namespace fast_1_3;


#define DEBUG

/*	Načtu cestu k síti
 *  Vytvořím prázdnou sít
 *  Přečtu sít pomoc GmshMeshReader
 *  Naplnim sít
 *  Vyšetřuji průniky
 * */

TEST(intersections, 1d_3d_fast){
	FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
	FilePath mesh_file("mesh/pokus5.msh", FilePath::input_file);
	Mesh mesh_krychle;

	ifstream ifs(string(mesh_file).c_str());

	mesh_krychle.read_gmsh_from_stream(ifs);

	InspectElements inspectelements( &mesh_krychle);

	vector<IntersectionLocal> il = inspectelements.getIntersections();

	xprintf(Msg,"PRUNIKU: %d \n", il.size());
	xprintf(Msg, "FAST INTERSECTION: Test is complete");
	inspectelements.print("aaa.txt");
}
