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
	FilePath mesh_file("mesh/tetra_abs.msh", FilePath::input_file); // krychle 1x1x1 param = 0.2; sít úseček param = 0.1
	Mesh mesh_krychle;

	ifstream ifs(string(mesh_file).c_str());

	mesh_krychle.read_gmsh_from_stream(ifs);

	//GmshMeshReader reader(mesh_file);

	//reader.read_mesh(&mesh_krychle);

	/*double *neco;
	double ok;

	*neco = 6.0;
	ok = *neco;

	xprintf(Msg, "ok:%f neco:%f \n", ok, *neco);*/


	InspectElements inspectelements( &mesh_krychle);

	vector<IntersectionLocal> il = inspectelements.getIntersections();

	xprintf(Msg,"PRUNIKU: %d \n", il.size());
	xprintf(Msg, "FAST INTERSECTION: Test is complete");

}




