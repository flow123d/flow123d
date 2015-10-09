/*
 * ProlongationLine.h
 *
 *  Created on: 2.12.2014
 *      Author: viktor
 */

#ifndef PROLONGATIONLINE_H_
#define PROLONGATIONLINE_H_

#include "plucker.h"

using namespace std;
namespace computeintersection{

/**
 * Simple class defines indices of elements for later processing of computing intersection 2D-3D
 * 
 * TODO:
 * - comment
 * - split prolongationline into Prolongation2D (neighbor is 2D) and Prolongation3D (neighbor is 3D)
 * - replace with struct, remove .cpp file
 */
class ProlongationLine {
private:

	// idx sousedniho elementu - podle typu je to 2D
	// elemnt nebo 3D element

	unsigned int elm_2D_idx;
	unsigned int elm_3D_idx;
	unsigned int dictionary_idx;

	int elm_2D_idx_old;
	int elm_3D_idx_old;


	// index of side -> takes C23->C13[:]->C12[side_idx]
	// every Plucker products
	// using ref simplex -> set to every new C12

	unsigned int side_idx_old;
	unsigned int side_idx;

	// for case side prolongation, there will be
	// always 9 Plucker products
	std::vector<double> pluckerProducts;

	// for case side prolongation, there will be
	// always 3 PC for 3D side and 3 PC for triangle
	std::vector<Plucker> pluckers;


	// for case triangle prolongation, there will be
	// only one PC for one edge of triangle
	// and 6 Plucker products

public:
	ProlongationLine(unsigned int element_2D,unsigned int element_3D, unsigned int dictionary, int element_2D_old = -1, int element_3D_old = -1);
	ProlongationLine();

	inline void setPluckerProducts(double pp, unsigned int index){
		pluckerProducts[index] = pp;
	};

	inline void setPlucker(Plucker p, unsigned int index){
		pluckers[index] = p;
	};

	inline double getPluckerProducts(unsigned int index) const{
		return pluckerProducts[index];
	};

	inline Plucker getPlucker(unsigned int index) const{
		return pluckers[index];
	};

	inline ~ProlongationLine(){};

	inline unsigned int getDictionaryIdx() const{
		return dictionary_idx;
	};

	inline unsigned int getElement2DIdx() const{
		return elm_2D_idx;
	};

	inline unsigned int getElement3DIdx() const{
		return elm_3D_idx;
	};

};

}
#endif /* PROLONGATIONLINE */
