/*
 * inspectelements.h
 *
 *  Created on: 11.3.2013
 *      Author: viktor
 *
 *      Uchováva informaci o všech elementech v sítí, jestli jsme je již prošli nebo ne
 *      Do konstruktoru přidá velikost sítě všech elementů
 *      metoda set_element nastaví konkrétní element, že jsme ho prošli
 *      metoda projety_element vrací informaci, jestli jsme už prošli daný element
 */

#include <vector>

class InspectElements {
private:
	vector<bool> projeti;

public:
	inline InspectElements(unsigned int size)
	: projeti(size, false)
	{
		//projeti.reserve(size);
		/*
		 * init:
		 * for(unsigned i = 0; i < size; i++){
		 * projeti[i] = false;
		 * }
		 *
		 * */
	}
	inline ~InspectElements(){}

	inline void set_element(unsigned int index){
		projeti[index] = true;
	}
	inline bool projety_element(unsigned int index){
		return projeti[index];
	}
};
