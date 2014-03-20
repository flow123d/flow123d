
#include <flow_gtest.hh>
//#define Flow123d_DEBUG
#include "system/sys_profiler.hh"

//#include "mesh/mesh.h"
//#include "mesh/msh_gmshreader.h"
//#include "mesh/bih_tree.hh"


#include "intersection/computeintersection.h"

using namespace std;
using namespace computeintersection;

TEST(intersections, all) {



	cout << "================== Testy objektů a potřebných metod ==================" << endl;
	cout << "========== Includované objekty: =============" << endl;
	arma::vec3 vectorA, vectorB, vectorC, vectorD;
	vectorA[0] = 1;vectorA[1] = 2;vectorA[2] = 3;
	vectorB[0] = 4;vectorB[1] = 15;vectorB[2] = 6;
	vectorC[0] = 27;vectorC[1] = 18;vectorC[2] = 39;
	vectorD[0] = 105;vectorD[1] = 111;vectorD[2] = 512;
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


	arma::vec3 vector;
	vector[0] = 1;
	vector[1] = 2;
	vector[2] = 3;
	arma::vec3 vector2;
	vector2[0] = 4;
	vector2[1] = 5;
	vector2[2] = 6;
	arma::vec3 vector3;
	vector3[0] = 7;
	vector3[1] = 8;
	vector3[2] = 9;
	arma::vec3 vector4;
	vector4[0] = 10;
	vector4[1] = 11;
	vector4[2] = 12;
	arma::vec3 pole_vectoru[] = {vector,vector2,vector3,vector4};
	arma::vec3 pole_vectoru2[] = {vectorB,vectorC,vectorD};


	Simplex<2> sim2(pole_vectoru2);
	Simplex<3> sim3(pole_vectoru);

	ComputeIntersection<Simplex<2>, Simplex<3> > novyCI(sim2, sim3);

	Plucker *muj_plucker = new Plucker(vector, vector2);
	Plucker ha(vector, vector3);
	Plucker da(vector, vector4);

	cout << "Testovaci Plucker[2]: ";
	ha.toString();

	novyCI.setPC_tetrahedron(&ha, 2);

	cout << "Zkouska1 Plucker[2]: "; (*novyCI.getPC_tetrahedron(2)).toString();

	novyCI.init();

	xprintf(Msg, "hoho");
	cout << "==================================================================" << endl;
}




