/*
 * unit.cc
 *
 *  Created on: Aug 30, 2012
 *      Author: jb
 */


#include <sstream>

#include "fields/unit.hh"


using namespace std;


Unit::Unit() {
	exponents_.resize(7);
	std::fill(exponents_.begin(), exponents_.end(), 0);
}

const std::vector<std::string> Unit::unit_symbols={"m","kg","s","A","K","mol","cd" };

Unit & Unit::m(int val) {
	exponents_[0] = val;
	return *this;
}

Unit & Unit::kg(int val) {
	exponents_[1] = val;
	return *this;
}

Unit & Unit::s(int val) {
	exponents_[2] = val;
	return *this;
}

Unit & Unit::A(int val) {
	exponents_[3] = val;
	return *this;
}

Unit & Unit::K(int val) {
	exponents_[4] = val;
	return *this;
}

Unit & Unit::mol(int val)  {
	exponents_[5] = val;
	return *this;
}

Unit & Unit::cd(int val) {
	exponents_[6] = val;
	return *this;
}

std::string Unit::print() {
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
