
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
	arma::vec3 pole_vectoru2[] = {vector,vector2,vector3};


	Simplex<2> sim2(pole_vectoru2);
	Simplex<3> sim3(pole_vectoru);

	sim3.toString();
	sim2.toString();

	ComputeIntersection<Simplex<2>, Simplex<3> > novyCI2;
	ComputeIntersection<Simplex<2>, Simplex<3> > novyCI(sim2, sim3);

	Plucker *muj_plucker = new Plucker(vector, vector2);
	Plucker ha(vector, vector3);
	Plucker da(vector, vector4);

	cout << "Testovaci Plucker[2]: " << ha[2] << endl;

	novyCI.setPC_tetrahedron(&ha, 2);

	Plucker a = *novyCI.p_coordinates_tetrahedron[2];
	cout << "Zkouska1 Plucker[2]: " << a[2] << endl;

	xprintf(Msg, "hoho");
	cout << "==================================================================" << endl;
}




