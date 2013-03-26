/*
 * Intersection0D.h
 *
 *  Created on: 11.3.2013
 *      Author: viktor
 */

class Intersection0D {
	unsigned int index_1D_element;
	unsigned int index_3D_element;
	arma::vec3 prusecik;
public:
	inline Intersection0D(unsigned int index_1D, unsigned int index_3D, arma::vec3 prus):index_1D_element(index_1D), index_3D_element(index_3D),prusecik(prus){};
		inline ~Intersection0D(){};


	inline unsigned int getIndex1DElement() const{
	    return index_1D_element;}

	inline void setIndex1DElement(unsigned int index1DElement){
	    index_1D_element = index1DElement;}

	inline unsigned int getIndex3DElement() const{
	    return index_3D_element;}

	inline void setIndex3DElement(unsigned int index3DElement){
	    index_3D_element = index3DElement;}

	inline arma::vec3 getPrusecik(){
		return prusecik;
	}
};
