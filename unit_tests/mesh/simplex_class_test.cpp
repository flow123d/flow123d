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

TEST(simplex, all) {

	//FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

	//Profiler::initialize();


	arma::vec3 Point0;Point0[0] = 1;Point0[1] = 2;Point0[2] = 3;
	arma::vec3 Point1;Point1[0] = 2;Point1[1] = 3;Point1[2] = 4;
	arma::vec3 Point2;Point2[0] = 3;Point2[1] = 4;Point2[2] = 1;
	arma::vec3 Point3;Point3[0] = 4;Point3[1] = 1;Point3[2] = 2;
	//arma::vec3 Point4;Point4[0] = 5;Point4[1] = 5;Point4[2] = 5;


	arma::vec3 *pp, *ll;
	arma::vec3 *pole_p[4];
	//pp = new arma::vec3;
	pp = &Point0;
	ll = &Point0;

	arma::vec3 *pole_pp[] = {&Point0,&Point1,&Point2,&Point3};
	arma::vec3 pole[] = {Point0,Point1,Point2,Point3};
	pole_p[0] = &Point0;
	pole_p[1] = &Point1;
	pole_p[2] = &Point2;
	pole_p[3] = &Point3;
	       //= {pp,pp,pp,ll};

	Point0[0] = 10;

	(*pp)[0] = 20;

	(pole_p[0][0])[0] = 30;



	xprintf(Msg, "pp:%f, ll:%f\n",(*pp)[0],(*ll)[0]);
	//*pp0 = Point0;

	//Simplex<3> sss(pole);
	Simplex<3> ss(pole_pp);


	//sss.toString();
	ss.to_string();

	//Simplex<3> a;
	//a[0][0][0].setPointCoordinates(Point0);
	//a[0][0][1].setPointCoordinates(Point1);
	//a[0][1][1].setPointCoordinates(Point2);
	//a[1][1][1].setPointCoordinates(Point3);

	//a.to_string();


	for(unsigned int i = 0; i < 100000;i++){
		//Simplex<3> sss(pole);
		Simplex<3> ss(pole_pp);
	}
	//Profiler::uninitialize();

	xprintf(Msg, "Test complete!\n");
}




