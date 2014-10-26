#include <flow_gtest.hh>
//#define Flow123d_DEBUG
#include "system/system.hh"
#include "system/sys_profiler.hh"
#include "system/file_path.hh"
#include <array>
//#include "mesh/mesh.h"
#include "mesh/msh_gmshreader.h"
//#include "mesh/bih_tree.hh"
//#include "mesh/ngh/include/point.h"
//#include "mesh/ngh/include/intersection.h"


//#include "intersection/computeintersection.h"
#include "intersection/inspectelements.h"

using namespace std;
using namespace computeintersection;


TEST(intersections, all) {



	/*cout << "================== Testy objektů a potřebných metod ==================" << endl;
	cout << "========== Includované objekty: =============" << endl;
	arma::vec3 vectorA, vectorB, vectorC, vectorD;
	vectorA[0] = 1;vectorA[1] = 2;vectorA[2] = 3;
	vectorB[0] = 0;vectorB[1] = 0;vectorB[2] = 0;
	vectorC[0] = 10;vectorC[1] = 0;vectorC[2] = 0;
	vectorD[0] = 5;vectorD[1] = 10;vectorD[2] = 0;
	cout << "=== Vector(" << vectorA[0] << "," << vectorA[1] << "," << vectorA[2] << ") ===" << endl;
	arma::vec3 pole_vec[] = {vectorA,vectorB,vectorC,vectorD};


	cout << "========== Vlastní objekty: =============" << endl;

	Plucker test_plucker(vectorC, vectorD);
	Plucker test_plucker2(vectorA, vectorB);
	cout << "====== Objekt Plucker =====" << endl;
	cout << "=== .toString():";
		test_plucker2.toString();
	cout << "=== .operator*(double): ";
		test_plucker2*5;
		test_plucker2.toString();
	cout << "=== .operator*(Plucker): " << test_plucker*test_plucker2 << endl;


	Simplex<3> simplex_3(pole_vec);
	cout << "====== Objekt Simplex =====" << endl;
	cout << "=== .toString(): " << endl;
		simplex_3[0][0].toString();
	cout << "=== .getAbscissa(unsigned int): " << endl;
	simplex_3.getAbscissa(2);
	simplex_3[0].getAbscissa(1);
	simplex_3[0][0].getAbscissa(1).toString();
	cout << "=========================================" << endl;
	cout << "=========================================" << endl;
	cout << "=========================================" << endl;

	IntersectionLocal il2;*/
	/* ctyrsten uvnitr celeho trojuhelniku
	arma::vec3 bodA; bodA[0] = -5; bodA[1] = 1; bodA[2] = -5;
	arma::vec3 bodB; bodB[0] = 10; bodB[1] = 2; bodB[2] = -4;
	arma::vec3 bodC; bodC[0] = 6; bodC[1] = 3; bodC[2] = 10;
	*/
	/* jeden vrchol trojuhelnika konci v ctyrstenu*/
	/*arma::vec3 bodA; bodA[0] = -5; bodA[1] = 1; bodA[2] = -5;
	arma::vec3 bodB; bodB[0] = 10; bodB[1] = 2; bodB[2] = -4;
	arma::vec3 bodC; bodC[0] = 5.5; bodC[1] = 3; bodC[2] = 0.5;*/

	// Trojúhelník je uvnitř čtyřstěnu
	/*arma::vec3 bodA; bodA[0] = 2; bodA[1] = 0.5; bodA[2] = 1;
	arma::vec3 bodB; bodB[0] = 6; bodB[1] = 0.6; bodB[2] = 0.5;
	arma::vec3 bodC; bodC[0] = 5; bodC[1] = 0.2; bodC[2] = 2;


	arma::vec3 Alfa; Alfa[0] = 0; Alfa[1] = 0; Alfa[2] = 0;
	arma::vec3 Beta; Beta[0] = 10; Beta[1] = 0; Beta[2] = 0;
	arma::vec3 Gama; Gama[0] = 5; Gama[1] = 0; Gama[2] = 5;
	arma::vec3 Delta; Delta[0] = 5; Delta[1] = 5; Delta[2] = 0;

	arma::vec3 vec_2[] = {bodA, bodB, bodC};
	arma::vec3 vec_3[] = {Alfa, Beta, Gama, Delta};
	Simplex<2> sim_2(vec_2);
	Simplex<3> sim_3(vec_3);
	cout << "=================" << endl;
	ComputeIntersection<Simplex<2>, Simplex<3>> CI23(sim_2, sim_3);
	CI23.init();
	CI23.compute(il2);

	for(unsigned int i = 0; i < il2.getIPsize();i++){
		IntersectionPoint<2,3> IP23 = il2.get_point(i);
		IP23.getLocalCoords1().print();
		sim_2.toString();
		arma::vec3 T_globalni = (IP23.getLocalCoords1())[0] * sim_2[0][0].getPointCoordinates()
							   +(IP23.getLocalCoords1())[1] * sim_2[0][1].getPointCoordinates()
							   +(IP23.getLocalCoords1())[2] * sim_2[1][1].getPointCoordinates();
		arma::vec3 C_globalni = (IP23.getLocalCoords2())[0] * sim_3[0][0][0].getPointCoordinates()
							   +(IP23.getLocalCoords2())[1] * sim_3[0][0][1].getPointCoordinates()
							   +(IP23.getLocalCoords2())[2] * sim_3[0][1][1].getPointCoordinates()
							   +(IP23.getLocalCoords2())[3] * sim_3[1][1][1].getPointCoordinates();
		cout << "musi byt stejne - tenhle:" << endl;
		T_globalni.print();
		cout << "a tento:" << endl;
		C_globalni.print();
	}*/


	cout << "===============" << endl;
	//FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
	FilePath mesh_file("/home/viktor/diplomka/flow123d/unit_tests/mesh/site/triangle_tetrahedron11.msh", FilePath::input_file);

	Profiler::initialize();


	Mesh mesh;

	GmshMeshReader reader(mesh_file);

	reader.read_mesh(&mesh);


	 cout << "Síť načtena!" << endl;
	cout << "Probíhá výpočet průniku" << endl;

	//Profiler::initialize();

	InspectElements ie(&mesh);
	ie.print(0);
	ie.print(1);
	double obsah = ie.polygonArea();
	xprintf(Msg,"Obsah polygonu: %f\n", obsah);
	//Profiler::instance()->output(0,cout);//MPI_COMM_WORLD,cout);
	Profiler::uninitialize();

	xprintf(Msg, "test nulových pluckerových souřadnic");





	xprintf(Msg, "Test complete!");
}




