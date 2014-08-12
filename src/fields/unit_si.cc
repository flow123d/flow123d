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

UnitSI & UnitSI::m(int val) {
	exponents_[0] = val;
	return *this;
}

UnitSI & UnitSI::kg(int val) {
	exponents_[1] = val;
	return *this;
}

UnitSI & UnitSI::s(int val) {
	exponents_[2] = val;
	return *this;
}

UnitSI & UnitSI::A(int val) {
	exponents_[3] = val;
	return *this;
}

UnitSI & UnitSI::K(int val) {
	exponents_[4] = val;
	return *this;
}

UnitSI & UnitSI::mol(int val)  {
	exponents_[5] = val;
	return *this;
}

UnitSI & UnitSI::cd(int val) {
	exponents_[6] = val;
	return *this;
}

std::string UnitSI::print() {
	std::stringstream output;
	for (unsigned int i=0; i<exponents_.size(); i++)
		if (exponents_[i]) {
			if (output.str().size()>0) {
				output << ".";
			}

			if (exponents_[i] == 1) {
				output << unit_symbols[i];
			} else if (exponents_[i] < 0) {
				output << unit_symbols[i] << "^(" << exponents_[i] << ")";
			} else {
				output << unit_symbols[i] << "^" << exponents_[i];
			}
		}

	return output.str();
}
