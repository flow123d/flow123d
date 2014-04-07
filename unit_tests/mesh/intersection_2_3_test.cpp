
#include <flow_gtest.hh>
//#define Flow123d_DEBUG
#include "system/sys_profiler.hh"

//#include "mesh/mesh.h"
//#include "mesh/msh_gmshreader.h"
//#include "mesh/bih_tree.hh"


#include "intersection/computeintersection.h"

using namespace std;
using namespace computeintersection;


void neco(Plucker &p_ref){
	//cout << "adresa vlozeneho parametru: " << p_ref << endl;
	cout << "adresa vlozeneho parametru: " << &p_ref << endl;

	Plucker *tak = &p_ref;
	cout << "adresa noveho pointru - na co koukam:" << tak << endl;
	cout << "adresa pointu: " << &tak << endl;
	tak = new Plucker();
	cout << "adresa noveho pointru - na co koukam:" << tak << endl;
	//Plucker ha = p_ref;
	//cout << "Adresa noveho objektu: " << &ha << endl;

}

TEST(intersections, all) {



	cout << "================== Testy objektů a potřebných metod ==================" << endl;
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


	arma::vec3 vector;
	vector[0] = 1;
	vector[1] = 1;
	vector[2] = 1;
	arma::vec3 vector2;
	vector2[0] = 9;
	vector2[1] = 1;
	vector2[2] = 2;
	arma::vec3 vector3;
	vector3[0] = 5;
	vector3[1] = 9;
	vector3[2] = 3;
	arma::vec3 vector4;
	vector4[0] = 6;
	vector4[1] = 5;
	vector4[2] = -2;
	arma::vec3 pole_vectoru[] = {vector,vector2,vector3,vector4};
	arma::vec3 pole_vectoru2[] = {vectorB,vectorC,vectorD};

	Simplex<2> sim2(pole_vectoru2);
	Simplex<3> sim3(pole_vectoru);

	ComputeIntersection<Simplex<2>, Simplex<3> > novyCI(sim2, sim3);

	//novyCI.setPC_tetrahedron(test_plucker, 1);
	//novyCI.setPC_triangle(test_plucker2,0);
	novyCI.init();
	novyCI.compute();

	arma::vec3 bodA; bodA[0] = 3; bodA[1] = 3; bodA[2] =4;
	arma::vec3 bodB; bodB[0] = 3; bodB[1] = 3; bodB[2] =-4;
	arma::vec3 Alfa; Alfa[0] = 1; Alfa[1] = 2; Alfa[2] =4;
	arma::vec3 Beta; Beta[0] = 5; Beta[1] = 2; Beta[2] =4;
	arma::vec3 Gama; Gama[0] = 1; Gama[1] = 6; Gama[2] =4;


	arma::vec3 vec_1[] = {bodA, bodB};
	arma::vec3 vec_2[] = {Alfa, Beta, Gama};
	Simplex<1> sim_1(vec_1);
	Simplex<2> sim_2(vec_2);
	//ComputeIntersection<Simplex<1>, Simplex<2>> CI12(sim_1, sim_2);
	//CI12.compute();
	//CI12.toStringPluckerCoordinates();

	//novyCI.clear_all();
		//novyCI.toStringPluckerCoordinatesTree();
		/*
		arma::vec3 vecAA; vecAA[0] = 0; vecAA[1]= 2;vecAA[2] = 0;
		arma::vec3 vecBB; vecBB[0] = 2; vecBB[1]= 3;vecBB[2] = 4;
		arma::vec3 vecCC; vecCC[0] = 5; vecCC[1]= 1;vecCC[2] = -1;
		arma::vec3 vecDD; vecDD[0] = 6; vecDD[1]= 6;vecDD[2] = 1;

		Plucker AB(vecAA, vecBB);
		Plucker CD(vecCC, vecDD);
		Plucker DC(vecDD, vecCC);

		double soucinABCD = AB*CD;
		double soucinABDC = AB*DC;

		AB.toString();
		CD.toString();
		DC.toString();
		cout << soucinABCD << " - " << soucinABDC << endl;
		cout << CD*AB << "-" << AB*CD << endl;
		cout << DC*AB << "-" << AB*DC << endl;


		cout << "test referenci: " << endl;
		AB.toString();
		cout << "adresa: " << &AB << endl;

		Plucker *kam_koukam = NULL;
		cout << "adresa kam koukam uvnitr pointru: " << kam_koukam << endl;
		cout << "primo adresa pointru:" << &kam_koukam << endl;
		neco(*kam_koukam);

*/

	cout << "==================================================================" << endl;
}




