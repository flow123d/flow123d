// bp_fris.cpp : Defines the entry point for the console application.
//
#include <stdio.h>
#include <iostream>
#include <time.h>
#include "intersections.h"
#include <vector>
using namespace std;

int main() {
	cout << "ZKOUSKA VSECH TRID" << endl;
	cout << "==================" << endl;
	cout << "BODY:" << endl;
	cout << "=====" << endl;
	Point<3> a(0, 0, 0);
	a.toString();
	Point<3> b(0, 0, 1);
	b.toString();
	Point<3> c(1, 1, 1);
	c.toString();
	Point<3> d(1, 1, 0);
	d.toString();
	Point<3> e(2, 5, 4);
	e.toString();

	Vector<3> aa(a, b);
	Vector<3> bb(d, e);

	Plucker pl(aa, bb);

	HyperPlane<1, 3> hp(Point<3>(1, 1, 0), Point<3>(1, 1.1, 3));
	HyperPlane<0, 3> hp0(b);

	Simplex<0, 3> s0(a);
	Point<3> pole[] = { a, b };
	Simplex<1, 3> s1(pole);
	Point<3> pole2[] =
			{ Point<3>(0, 0, 0), Point<3>(3, 0, 0), Point<3>(0, 3, 0) };
	Simplex<2, 3> s2(pole2);
	Point<3> pole3[] = { a, b, c, d };
	Simplex<3, 3> s3(pole3);

	Point<3> prusecik;
	prusecik = Intersections(hp, s2);
	cout << "prusecik lokálně:";
	prusecik.toString();
	prusecik = Globalni_souradnice(prusecik, s2);
	cout << "globální: ";
	prusecik.toString();

	/*Pnn = Intersections(hp1, moje[0]);
	 Pnn->toString();
	 Point<3> *pplo;
	 pplo = Intersections(hp1, moje[0]);
	 if(pplo == NULL){cout << "yes";}
	 else{cout << "no";}
	 pplo = &ee;
	 if(pplo == NULL){cout << "yes";}
	 else{cout << "no";}
	 */
	//moje.toString();
	//moje.Simplex()
	//seru = Simplex<1,3>(koko);
	//Simplex<1,3> muj = new Simplex<1,3>(koko);

	/*
	 srand(time(NULL));
	 vector<Bod<3>> body_site;												//Vytvo�en� pole bod�
	 body_site.reserve(10000);												//alokovano misto
	 for(int i = 0;i< 10000;i++){
	 body_site.push_back(Bod<3>(rand()%10000,rand()%10000,rand()%10000));	//vytvoreni nahodneho bodu
	 }

	 vector< Simplex<3,3> > simplices_3;									//vytvoreni pole simplexu3
	 simplices_3.reserve(10000);												//alokovano pole
	 for(int i=0; i< 10000; i++) {
	 Bod<3> ole[] = {body_site[rand()%10000],body_site[rand()%10000],body_site[rand()%10000],body_site[rand()%10000]};
	 simplices_3.push_back(Simplex<3,3>(ole));
	 }


	 double soucet = 0;
	 for(int k = 0; k < 10; k++){
	 HyperPlane<1,3> mojeprimka(Bod<3>(rand()%10000,rand()%10000,rand()%10000), Bod<3>(rand()%10000,rand()%10000,rand()%10000));		// libovolna primka
	 Bod<3> prus(NULL,NULL,NULL), prus2(NULL,NULL,NULL);

	 double ted = timeGetTime();
	 for(int j = 0; j <10000;j++){
	 Intersections(mojeprimka, simplices_3[j], prus, prus2);
	 prus = Bod<3>(NULL,NULL,NULL);prus2 = Bod<3>(NULL,NULL,NULL);
	 }
	 double potom = timeGetTime() - ted;
	 cout << "Zak: " << potom << endl;

	 double ted2 = timeGetTime();
	 for(int j = 0; j < 10000; j++){
	 IntersectionsOp(mojeprimka, simplices_3[j], prus, prus2);
	 prus = Bod<3>(NULL,NULL,NULL);prus2 = Bod<3>(NULL,NULL,NULL);
	 }
	 double potom2 = timeGetTime() - ted2;
	 cout << "Opt: " << potom2 << endl;

	 soucet += (100*potom/potom2 - 100);
	 cout << "Je rychlejsi o: " << 100*potom/potom2 - 100 << "%" << endl;
	 }
	 cout << "pr�m�rn� je optimalizace lep�� o: " << soucet/10 << "%" << endl;
	 */

	return 0;
}

