/*
 * ProlongationPoint1D.h
 *
 *  Created on: 20.3.2013
 *      Author: viktor
 */

class PPBase;

class ProlongationPoint1D : public PPBase{
	unsigned int vychozi3Delement;
public:
	inline ProlongationPoint1D(unsigned int index):vychozi3Delement(index){};
	inline ~ProlongationPoint1D(){};
	inline void VypocitejPrvni3D(){
		/* Díky informaci o přechozím 1D elementu a další stěně vypočítá navazující 1D element
		 * S ním a s výchozím 3D elementem spočítá počáteční 3D prolongation point
		 * označí 1D i 3D jako projitý
		 *
		 * /

	}
};
