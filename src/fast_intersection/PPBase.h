/*
 * PPBase.h
 *
 *  Created on: 20.3.2013
 *      Author: viktor
 */

class Intersection0D;

class PPBase {
	Element elm; //nebo index elementu
	Intersection0D prusecik; // U 1D elementu to bude vlastně počáteční a koncový bod
	Side stena;
public:
	inline PPBase(){};
	inline ~PPBase(){};
};

