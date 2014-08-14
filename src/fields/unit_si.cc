/*
 * unit_si.cc
 *
 *  Created on: Aug 30, 2012
 *      Author: jb
 */


#include <sstream>

#include "fields/unit_si.hh"
#include "system/xio.h"


using namespace std;


UnitSI::UnitSI() {
	exponents_.resize(7);
	std::fill(exponents_.begin(), exponents_.end(), 0);
	undef_ = true;
}

const UnitSI UnitSI::N = UnitSI().m().kg().s(-2);;
const UnitSI UnitSI::J = UnitSI().m(2).kg().s(-2);;
const UnitSI UnitSI::W = UnitSI().m(2).kg().s(-3);;
const UnitSI UnitSI::Pa = UnitSI().m(-1).kg().s(-2);;

const std::vector<std::string> UnitSI::unit_symbols={"m","kg","s","A","K","mol","cd" };

UnitSI & UnitSI::m(int exp) {
	exponents_[0] = exp;
	undef_ = false;
	return *this;
}

UnitSI & UnitSI::kg(int exp) {
	exponents_[1] = exp;
	undef_ = false;
	return *this;
}

UnitSI & UnitSI::s(int exp) {
	exponents_[2] = exp;
	undef_ = false;
	return *this;
}

UnitSI & UnitSI::A(int exp) {
	exponents_[3] = exp;
	undef_ = false;
	return *this;
}

UnitSI & UnitSI::K(int exp) {
	exponents_[4] = exp;
	undef_ = false;
	return *this;
}

UnitSI & UnitSI::mol(int exp)  {
	exponents_[5] = exp;
	undef_ = false;
	return *this;
}

UnitSI & UnitSI::cd(int exp) {
	exponents_[6] = exp;
	undef_ = false;
	return *this;
}

std::string UnitSI::print() const {
	ASSERT(is_def(), "UnitSI object must be defined!");
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

void UnitSI::undef(bool val) {
	undef_ = val;
}

bool UnitSI::is_def() const {
	return !undef_;
}


UnitSI operator *(const UnitSI &a, const UnitSI &b) {
	UnitSI tmp;

	if (a.is_def() && b.is_def()) {
		tmp.undef_ = false;
		for (unsigned int i=0; i<7; i++) {
			tmp.exponents_[i] = a.exponents_[i] + b.exponents_[i];
		}
	}

	return tmp;
}
