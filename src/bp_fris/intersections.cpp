#include "intersections.h"

SPoint<3> Intersections(HyperPlane<1,3> a, Simplex<2,3> b){
	double c,d,e;
	//Point<3> lokalni_souradnice();
	SPoint<3> ba = b[0][0].getPoint();			// Bod A
	SPoint<3> bb = b[0][1].getPoint();			// Bod B
	SPoint<3> bc = b[1][1].getPoint();			// Bod C
	Plucker pluk(Vector<3>(ba,bb),Vector<3>(ba,bb)*(Vector<3>)ba);  // sm�r AB
	Plucker pluk2(Vector<3>(bb,bc),Vector<3>(bb,bc)*(Vector<3>)bb);	// BC
	Plucker pluk3(Vector<3>(bc,ba),Vector<3>(bc,ba)*(Vector<3>)bc);	// CA

	c = a.getPlucker()*pluk;
	d = a.getPlucker()*pluk2;
	e = a.getPlucker()*pluk3;
	if((c > 0 && d > 0 && e > 0) || (c < 0 && d < 0 && e < 0)){

		SPoint<3> vysledek(d/(c+d+e),e/(c+d+e),c/(c+d+e));

	return vysledek;	// loukalni souradnice u0, u1, u2
	}
	else{
		/*if(c == 0){return Intersections(a,b[0]);}
		if(d == 0){return Intersections(a,b[1]);}
		if(e == 0){return Intersections(a,b[2]);}*/
	return SPoint<3>(0,0,0); // trojuhelnik se neprotne
	}
};


SPoint<3> Intersections(HyperPlane<1,3> a, Simplex<1,3> b){
	/* 
	Pr�se��k p��mka a �se�ka, parametrick� vyj�d�en� obou, �e�en� t�� rovnic o 2 nezn�m�ch
	kdy� jsou v�ety 't' stejn�, prot�naj� se, jinak jsou rovnob�n�, kdy� vektor je n�sobkem druh�ho
	jsou toto�n�... kdyz je t mezi 0 a 1, bod lezi na simplexu
	*/
	double a0, a1, a2, c0, c1, c2, x, y;
	Vector<3> u(b[0].getPoint(), b[1].getPoint());
	Vector<3> v(a.getPointA(), a.getPointB());
	double cislo[3];
	cislo[0] = b[0].getPoint()[0] - a.getPointA()[0];
	cislo[1] = b[0].getPoint()[1] - a.getPointA()[1];
	cislo[2] = b[0].getPoint()[2] - a.getPointA()[2];

	a0 = b[0].getPoint()[0];
	a1 = b[0].getPoint()[1];
	a2 = b[0].getPoint()[2];
	
	c0 = a.getPointA()[0];
	c1 = a.getPointA()[1];
	c2 = a.getPointA()[2];
	// totoznost nebo rovnobeznost:
	if(u[0] == v[0] && u[1] == v[1] && u[2] == v[2]){ // doplnit ze vektor u muze byt nasobkem vektoru v
		if( (cislo[0]/u[0]) == (cislo[1]/u[1]) == (cislo[2]/u[2])){ cout << "jsou totozne" << endl;} // !! o�et�it d�len� nulou
		else{cout << "jsou rovnobezne" <<endl;} 
	// �e�en� soustavy rovnic:
	} else if (u[0] != 0 && u[1] != 0 && u[2] != 0 && v[0] != 0 && v[1] != 0 && v[2] != 0){
		x = (v[0]*(c1-a1) + a0 - c0)/(v[0]*u[1] - u[0]);
		y = (a0 + u[0]*x - c0)/v[0];			// v�echny koeficienty nenulove, vypocet x + zkouska
		if((c2 - a2 + v[2]*y - u[2]*x) == 0 && (x >= 0 && x <= 1)){cout << "koeficient = " << x <<endl;}
		else{cout << "neprotina" <<endl;}
	} else if(v[0] == 0 && u[0] != 0){x = -cislo[0]/u[0];}
	else if(v[1] == 0 && u[1] != 0){x = -cislo[1]/u[1];}   // jeli jeden koef u Y nula, m�me X
	else if(v[2] == 0 && u[2] != 0){x = -cislo[2]/u[2];}
	else{
		if(v[0] != 0 && u[0] == 0){y = cislo[0]/v[0];}
		else if(v[1] != 0 && u[1] == 0){y = cislo[1]/v[1];} // jeli jeden koef u X nula, m�me Y
		else if(v[2] != 0 && u[2] == 0){y = cislo[2]/v[2];}
		// nemame ani X ani Y.. soustava 2 rovnic:
		else if(v[0] == 0 && u[0] == 0){x = (v[2]*cislo[1] - v[1]*cislo[2])/(v[1]*u[2] - v[2]*u[1]);}
		else if(v[1] == 0 && u[1] == 0){x = (v[2]*cislo[0] - v[0]*cislo[2])/(v[0]*u[2] - v[2]*u[0]);}
		else{x = (v[1]*cislo[0] - v[0]*cislo[1])/(v[0]*u[1] - v[1]*u[0]);}

		if(v[0] != 0 && u[0] != 0){x = (-cislo[0] + v[0]*y)/u[0];}
		else if(v[1] != 0 && u[1] != 0){x = (-cislo[1] + v[1]*y)/u[1];}
		else{x = (-cislo[2] + v[2]*y)/u[2];}
	}
	cout << "X: " << x <<endl;
	if(x >=0 && x <= 1){cout <<"primka protina simplex" << endl;}
	else{cout <<"primka neprotina simplex" << endl;}
	/*�e�en� soustavy 3 rovnic o 2 nezn�m�ch
	A0 + u0 . X = C0 + v0 . Y
	A1 + u1 . X = C1 + v1 . Y
	A2 + u2 . X = C2 + v2 . Y
	*/
return SPoint<3>(0,0,0);
};

void Intersections(HyperPlane<1,3> a, Simplex<3,3> b, SPoint<3> &prus1, SPoint<3> &prus2){
	for(int i =0; i < 4;i++){
		if(prus1[0] == 0){prus1 = Intersections(a,b[i]);}
		else if(prus2[0] == 0){prus2 = Intersections(a,b[i]);}
		else{break;} 
	}
};

bool Intersections(HyperPlane<1,3> a, Simplex<0,3> b){
return true;
};

double soucinAB, soucinBC, soucinCA, soucinBD, soucinDA, soucinCD;
SPoint<3>* IntersectionsOpp(HyperPlane<1,3> hp, Simplex<2,3> sm, int stena){
	double c,d,e;
	if(stena == 0){
		//AB = Plucker(Vektor<3>(sm[0][0].getBod(),sm[0][1].getBod()), Vektor<3>(sm[0][0].getBod(),sm[0][1].getBod())*(Vektor<3>)sm[0][0].getBod());
		//BC = Plucker(Vektor<3>(sm[0][1].getBod(),sm[1][1].getBod()), Vektor<3>(sm[0][1].getBod(),sm[1][1].getBod())*(Vektor<3>)sm[0][1].getBod());
		//CA = Plucker(Vektor<3>(sm[1][1].getBod(),sm[0][0].getBod()), Vektor<3>(sm[1][1].getBod(),sm[0][0].getBod())*(Vektor<3>)sm[1][1].getBod());
		c = soucinAB = hp.getPlucker()*Plucker(Vector<3>(sm[0][0].getPoint(),sm[0][1].getPoint()), Vector<3>(sm[0][0].getPoint(),sm[0][1].getPoint())*(Vector<3>)sm[0][0].getPoint());
		d = soucinBC = hp.getPlucker()*Plucker(Vector<3>(sm[0][1].getPoint(),sm[1][1].getPoint()), Vector<3>(sm[0][1].getPoint(),sm[1][1].getPoint())*(Vector<3>)sm[0][1].getPoint());
		e = soucinCA = hp.getPlucker()*Plucker(Vector<3>(sm[1][1].getPoint(),sm[0][0].getPoint()), Vector<3>(sm[1][1].getPoint(),sm[0][0].getPoint())*(Vector<3>)sm[1][1].getPoint());

	}
	else if(stena == 1){
		//AB;
		//BD = Plucker(Vektor<3>(sm[0][1].getBod(),sm[1][1].getBod()), Vektor<3>(sm[0][1].getBod(),sm[1][1].getBod())*(Vektor<3>)sm[0][1].getBod());
		//DA = Plucker(Vektor<3>(sm[1][1].getBod(),sm[0][0].getBod()), Vektor<3>(sm[1][1].getBod(),sm[0][0].getBod())*(Vektor<3>)sm[1][1].getBod());
		c = soucinAB;
		d = soucinBD = hp.getPlucker()*Plucker(Vector<3>(sm[0][1].getPoint(),sm[1][1].getPoint()), Vector<3>(sm[0][1].getPoint(),sm[1][1].getPoint())*(Vector<3>)sm[0][1].getPoint());
		e = soucinDA = hp.getPlucker()*Plucker(Vector<3>(sm[1][1].getPoint(),sm[0][0].getPoint()), Vector<3>(sm[1][1].getPoint(),sm[0][0].getPoint())*(Vector<3>)sm[1][1].getPoint());
	}
	else if(stena == 2){
		//CA;
		//CD = Plucker(Vektor<3>(sm[0][1].getBod(),sm[1][1].getBod()), Vektor<3>(sm[0][1].getBod(),sm[1][1].getBod())*(Vektor<3>)sm[0][1].getBod());
		//DA;
		c = -soucinCA;
		d = soucinCD = hp.getPlucker()*Plucker(Vector<3>(sm[0][1].getPoint(),sm[1][1].getPoint()), Vector<3>(sm[0][1].getPoint(),sm[1][1].getPoint())*(Vector<3>)sm[0][1].getPoint());
		e = soucinDA;
	}
	else{
	//BC;
	//CD;
	//BD; 
	c = soucinBC;
	d = soucinCD;
	e = -soucinBD;
	}
	
	//Bod<3> lokalni_souradnice();
	
	if((c > 0 && d > 0 && e > 0) || (c < 0 && d < 0 && e < 0)){
		SPoint<3> sp(d/(c+d+e),e/(c+d+e),c/(c+d+e));
		return &sp;	// loukalni souradnice u0, u1, u2
	}
	else{
	
		return NULL;//SPoint<3>(0,0,0);
				}
};

double productAB, productBC, productCA, productBD, productDA, productCD;
bool IntersectionsOp_1D_2D(HyperPlane<1,3> &hp, Simplex<2,3> sm, int stena, std::vector<double> &coords_3D, double &local_abscissa, bool &orientace){
	double c,d,e;
	if(stena == 0){
		c = productAB = hp.getPlucker()*Plucker(Vector<3>(sm[0][0].getPoint(),sm[0][1].getPoint()), Vector<3>(sm[0][0].getPoint(),sm[0][1].getPoint())*(Vector<3>)sm[0][0].getPoint());
		d = productBC = hp.getPlucker()*Plucker(Vector<3>(sm[0][1].getPoint(),sm[1][1].getPoint()), Vector<3>(sm[0][1].getPoint(),sm[1][1].getPoint())*(Vector<3>)sm[0][1].getPoint());
		e = productCA = hp.getPlucker()*Plucker(Vector<3>(sm[1][1].getPoint(),sm[0][0].getPoint()), Vector<3>(sm[1][1].getPoint(),sm[0][0].getPoint())*(Vector<3>)sm[1][1].getPoint());

	}
	else if(stena == 1){
		c = -productAB;
		productBD = hp.getPlucker()*Plucker(Vector<3>(sm[0][1].getPoint(),sm[1][1].getPoint()), Vector<3>(sm[0][1].getPoint(),sm[1][1].getPoint())*(Vector<3>)sm[0][1].getPoint());
		d = -productBD;
		productDA = hp.getPlucker()*Plucker(Vector<3>(sm[1][1].getPoint(),sm[0][0].getPoint()), Vector<3>(sm[1][1].getPoint(),sm[0][0].getPoint())*(Vector<3>)sm[1][1].getPoint());
		e = -productDA;
	}
	else if(stena == 2){
		c = -productCA;
		d = productCD = hp.getPlucker()*Plucker(Vector<3>(sm[0][1].getPoint(),sm[1][1].getPoint()), Vector<3>(sm[0][1].getPoint(),sm[1][1].getPoint())*(Vector<3>)sm[0][1].getPoint());
		e = productDA;
	}
	else{
	c = -productBC;
	d = -productCD;
	e = productBD;
	}

	/* Při výpočtu jsou P. souřadnice zorientovány, aby vycházely kladné hodnoty, pokud úsečka vstupuje do 4 stěnu
	 * a záporné pokud vystupuje
	 * */

	if((c > 0 && d > 0 && e > 0) || (c < 0 && d < 0 && e < 0)){
		// c = w0; d = w1; e = w2
		// lokální alfa = w2/soucet; lokální beta = w1/soucet; => lokální souřadnice na stěně
		double alfa = e/(c+d+e);
		double beta = c/(c+d+e);
		// lokální souřadnice na přímce T
		// T = localAbscissa= - A(i) + ( 1 - alfa - beta ) * V0(i) + alfa * V1(i) + beta * V2 (i) / U(i)
		// i = max z U(i)
		Vector<3> vec(hp.getPointA(),hp.getPointB());
		int i = 0;
		double max = vec[0];
		if(vec[1] > max) max = vec[1]; i = 1;
		if(vec[2] > max) max = vec[2]; i = 2;

		local_abscissa = (-hp.getPointA()[i] + (1 - alfa - beta)*sm[0][0].getPoint()[i] +
				alfa*sm[0][1].getPoint()[i] + beta*sm[1][1].getPoint()[i])/max;

		coords_3D.push_back(alfa); coords_3D.push_back(beta);

		if(c*d*e > 0){
			orientace = true; // orientace je v pořádku -> úsečka vchází do 4stěnu
		}
		else{
			orientace = false; // úsečka vychází
		}


		return true;
	}
	else{
		return false;
    }
};

void IntersectionsOp(HyperPlane<1,3> hp, Simplex<3,3> sm, SPoint<3> *prusecik1, SPoint<3> *prusecik2){
	for(int i = 0;i<4;i++){
		if(prusecik1 == NULL){prusecik1 = IntersectionsOpp(hp, sm[i], i);}
		else if(prusecik2 == NULL){prusecik2 = IntersectionsOpp(hp, sm[i], i);}
		else{break;}
	}
};

SPoint<3> Globalni_souradnice(SPoint<3> lokalni_souradnice, Simplex<2,3> trojuhelnik){

	SPoint<3> ba = trojuhelnik[0][0].getPoint();	// Bod A
	SPoint<3> bb = trojuhelnik[0][1].getPoint();			// Bod B
	SPoint<3> bc = trojuhelnik[1][1].getPoint();			// Bod C

	/* V�po�et glob�ln�ch sou�adnic X = u0 * V0 + u1 * V1 + u2 * V2 */
	return SPoint<3>((lokalni_souradnice[0]*ba[0] + lokalni_souradnice[1]*bb[0] + lokalni_souradnice[2]*bc[0]),
		(lokalni_souradnice[0]*ba[1] + lokalni_souradnice[1]*bb[1] + lokalni_souradnice[2]*bc[1]),
		(lokalni_souradnice[0]*ba[2] + lokalni_souradnice[1]*bb[2] + lokalni_souradnice[2]*bc[2]));
};
