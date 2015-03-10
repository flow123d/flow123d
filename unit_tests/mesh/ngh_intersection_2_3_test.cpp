#define TEST_USE_MPI
#include <flow_gtest_mpi.hh>
//#define Flow123d_DEBUG
#include "system/system.hh"
#include "system/sys_profiler.hh"
#include "system/file_path.hh"
#include <array>
#include "mesh/mesh.h"
#include "mesh/msh_gmshreader.h"
#include "mesh/bih_tree.hh"
#include "mesh/ngh/include/point.h"
#include "mesh/ngh/include/intersection.h"
#include "fields/field_interpolated_p0.hh"

using namespace std;


TEST(ngh_intersection_2_3, all) {

	unsigned int elementLimit = 20;
	cout << "===============" << endl;
	FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
	FilePath mesh_file("mesh/site/megasit10.msh", FilePath::input_file);

	Profiler::initialize();

	Mesh mesh;
	ifstream ifs(string(mesh_file).c_str());
	mesh.read_gmsh_from_stream(ifs);

	cout << "Síť načtena!" << endl;
	cout << "Probíhá výpočet průniku" << endl;

	{ START_TIMER("Vypocet pruniku");
	BIHTree bt(&mesh, elementLimit);

	double obsah = 0;

	FOR_ELEMENTS(&mesh, elm) {

		 if (elm->dim() == 2) {
			TTriangle tt;
			tt.SetPoints(TPoint(elm->node[0]->point()(0), elm->node[0]->point()(1), elm->node[0]->point()(2)),
									 TPoint(elm->node[1]->point()(0), elm->node[1]->point()(1), elm->node[1]->point()(2)),
									 TPoint(elm->node[2]->point()(0), elm->node[2]->point()(1), elm->node[2]->point()(2)) );
			std::vector<unsigned int> searchedElements;
			BoundingBox elementBoundingBox = tt.get_bounding_box();
			bt.find_bounding_box(elementBoundingBox, searchedElements);
			TTetrahedron th;
			TIntersectionType iType;
			double measure;

			{ START_TIMER("Hlavni vypocet");
			for (std::vector<unsigned int>::iterator it = searchedElements.begin(); it!=searchedElements.end(); it++)
				{
					int idx = *it;
					ElementFullIter ele = mesh.element( idx );
					if (ele->dim() == 3) {
						th.SetPoints(TPoint(ele->node[0]->point()(0), ele->node[0]->point()(1), ele->node[0]->point()(2)),
								TPoint(ele->node[1]->point()(0), ele->node[1]->point()(1), ele->node[1]->point()(2)),
								TPoint(ele->node[2]->point()(0), ele->node[2]->point()(1), ele->node[2]->point()(2)),
								TPoint(ele->node[3]->point()(0), ele->node[3]->point()(1), ele->node[3]->point()(2)) );

						//FieldInterpolatedP0<3,FieldValue<3>::Scalar>::createTetrahedron(ele, tt);
						GetIntersection(tt, th, iType, measure);
						obsah += measure;
						//if (iType == line) {rintf(Msg, "%d %d \n",elm.id(),ele.id()); }l
					}
				}
			END_TIMER("Hlavni vypocet");}
		 }
   }



	xprintf(Msg,"Obsah polygonu: %f\n", obsah);
	END_TIMER("Vypocet pruniku");}

	Profiler::instance()->output(MPI_COMM_WORLD, cout);
	Profiler::uninitialize();
	xprintf(Msg, "Test complete!");
}




