/*
 * unit_si.cc
 *
 *  Created on: Aug 30, 2012
 *      Author: jb
 */


#include <sstream>

#include "fields/unit_si.hh"


using namespace std;


UnitSI::UnitSI() {
	exponents_.resize(7);
	std::fill(exponents_.begin(), exponents_.end(), 0);
}

const std::vector<std::string> UnitSI::unit_symbols={"m","kg","s","A","K","mol","cd" };

UnitSI & UnitSI::m(int exp) {
	exponents_[0] = exp;
	return *this;
}

UnitSI & UnitSI::kg(int exp) {
	exponents_[1] = exp;
	return *this;
}

UnitSI & UnitSI::s(int exp) {
	exponents_[2] = exp;
	return *this;
}

UnitSI & UnitSI::A(int exp) {
	exponents_[3] = exp;
	return *this;
}

UnitSI & UnitSI::K(int exp) {
	exponents_[4] = exp;
	return *this;
}

UnitSI & UnitSI::mol(int exp)  {
	exponents_[5] = exp;
	return *this;
}

UnitSI & UnitSI::cd(int exp) {
	exponents_[6] = exp;
	return *this;
}

std::string UnitSI::print() {
	std::stringstream output;
	output << "$[";
	for (unsigned int i=0; i<exponents_.size(); i++)
		if (exponents_[i]) {
			if (exponents_[i] == 1) {
				output << unit_symbols[i];
			} else {
				output << unit_symbols[i] << "^{" << exponents_[i] << "}";
			}
		}

	if (output.str().size()==2) { //dimensionless quantity, contains only "$["
		output << "-";
	}
	output << "]$";

	return output.str();
}
