#define TEST_USE_MPI
#include <flow_gtest_mpi.hh>
//#define Flow123d_DEBUG
#include "system/system.hh"
#include "system/sys_profiler.hh"
#include "system/file_path.hh"
#include "intersection/inspectelements.h"

using namespace std;
using namespace computeintersection;


TEST(intersections, all) {

	cout << "===============" << endl;
	FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
	FilePath mesh_file("mesh/site/triangle_tetrahedron.msh", FilePath::input_file);
//     FilePath mesh_file("mesh/site/notfunctional/triangle_tetrahedron2.msh", FilePath::input_file);

	Profiler::initialize();

    unsigned int permutations[6][3] = {{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
    for(unsigned int p=0; p<6; p++)
    {
        Mesh mesh;
        // read mesh with gmshreader
        GmshMeshReader reader(mesh_file);
        reader.read_mesh(&mesh);
        
        // permute nodes:
        FOR_ELEMENTS(&mesh,ele)
        {
            if(ele->dim() == 2)
            {
                Node* tmp[3];
                for(unsigned int i=0; i<ele->n_nodes(); i++)
                {
                    tmp[i] = ele->node[permutations[p][i]];
                }
                for(unsigned int i=0; i<ele->n_nodes(); i++)
                {
                    ele->node[i] = tmp[i];
                    ele->node[i]->point().print(cout);
                }
            }
        }
        mesh.setup_topology();
        
        
        cout << "Síť načtena!" << endl;
        cout << "Probíhá výpočet průniku" << endl;

        InspectElements ie(&mesh);
        { START_TIMER("Vypocet pruniku");

        //ie.ComputeIntersections23();
        //ie.compute_intersections<2,4>();
        ie.compute_intersections<2,3>();
        END_TIMER("Vypocet pruniku");}
        //ie.print(0);
        //ie.print(1);
        ie.print_mesh_to_file("pp");

        double obsah = ie.polygonArea();
        xprintf(Msg,"Obsah polygonu: %f\n", obsah);
    }
    
	Profiler::instance()->output(MPI_COMM_WORLD, cout);
	Profiler::uninitialize();

	//count << pa.c_str() << endl;
	xprintf(Msg, "Test complete! ");
}




