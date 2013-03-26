/*
 * Intersection1D.h
 *
 *  Created on: 20.3.2013
 *      Author: viktor
 */


class Intersection1D {
	unsigned int prusecik1;
	unsigned int prusecik2;
public:
	inline Intersection1D(unsigned int prus1, unsigned int prus2):prusecik1(prus1), prusecik2(prus2){};
	inline ~Intersection1D(){};

	inline unsigned int getprusecik1(){
		return prusecik1;};
	inline unsigned int getprusecik2(){
		return prusecik2;};
};
