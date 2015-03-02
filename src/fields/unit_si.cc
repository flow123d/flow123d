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
	exponents_.resize(UnitSI::n_base_units);
	std::fill(exponents_.begin(), exponents_.end(), 0);
	undef_ = true;
}

UnitSI & UnitSI::N() {
	static UnitSI unit = UnitSI().m().kg().s(-2);
	return unit;
}

UnitSI & UnitSI::J() {
	static UnitSI unit = UnitSI().m(2).kg().s(-2);
	return unit;
}

UnitSI & UnitSI::W() {
	static UnitSI unit = UnitSI().m(2).kg().s(-3);
	return unit;
}

UnitSI & UnitSI::Pa() {
	static UnitSI unit = UnitSI().m(-1).kg().s(-2);
	return unit;
}

UnitSI & UnitSI::dimensionless() {
	static UnitSI unit = UnitSI().m(0);
	return unit;
}

UnitSI & UnitSI::m(int exp) {
	exponents_[UnitSI::order_m] = exp;
	undef_ = false;
	return *this;
}

UnitSI & UnitSI::kg(int exp) {
	exponents_[UnitSI::order_kg] = exp;
	undef_ = false;
	return *this;
}

UnitSI & UnitSI::s(int exp) {
	exponents_[UnitSI::order_s] = exp;
	undef_ = false;
	return *this;
}

UnitSI & UnitSI::A(int exp) {
	exponents_[UnitSI::order_A] = exp;
	undef_ = false;
	return *this;
}

UnitSI & UnitSI::K(int exp) {
	exponents_[UnitSI::order_K] = exp;
	undef_ = false;
	return *this;
}

UnitSI & UnitSI::mol(int exp)  {
	exponents_[UnitSI::order_mol] = exp;
	undef_ = false;
	return *this;
}

UnitSI & UnitSI::cd(int exp) {
	exponents_[UnitSI::order_cd] = exp;
	undef_ = false;
	return *this;
}

UnitSI & UnitSI::md(int exp) {
	exponents_[UnitSI::order_md] = exp;
	undef_ = false;
	return *this;
}






std::string UnitSI::format_latex() const {
    OutputFormat form;
    form.exp_open="^{";
    form.exp_close="}";
    form.delimiter="";
    return format(form);
}



std::string UnitSI::format_text() const {
    OutputFormat form;
    form.exp_open="(";
    form.exp_close=")";
    form.delimiter=".";
    return format(form);
}



std::string UnitSI::format(OutputFormat form) const {
	ASSERT(is_def(), "UnitSI object must be defined!");

    // Symbols for base SI units.
    std::vector<std::string> unit_symbols={"m","d","kg","s","A","K","mol","cd" };

	std::stringstream output;

	// format of meter (can be m^{n} or m^{n-d})
	if (exponents_[ UnitSI::order_m ] || exponents_[ UnitSI::order_md ]) {
		output << unit_symbols[ UnitSI::order_m ];
		if (exponents_[ UnitSI::order_m ]!=1 || exponents_[ UnitSI::order_md ]) {
			output << form.exp_open;
			if (exponents_[ UnitSI::order_m ]) {
				output << exponents_[ UnitSI::order_m ];
				if (exponents_[ UnitSI::order_md ]>0) output << "+";
			}
			if (exponents_[ UnitSI::order_md ]) {
				if (exponents_[ UnitSI::order_md ]==-1) output << "-";
				else if (exponents_[ UnitSI::order_md ]!=1) output << exponents_[ UnitSI::order_md ];
				output << unit_symbols[ UnitSI::order_md ];
			}
			output << form.exp_close;
		}
	}

	// format of other units
	for (unsigned int i=2; i<UnitSI::n_base_units; i++)
		if (exponents_[i]) {
		    if (output.str().size() > 0) output << form.delimiter;
		    output << unit_symbols[i];
			if (exponents_[i] != 1) output <<  form.exp_open << exponents_[i] << form.exp_close;
		}

	if (output.str().size()==0) { //dimensionless quantity, contains only "$["
		output << "-";
	}

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

UnitSI operator /(const UnitSI &a, const UnitSI &b) {
	UnitSI tmp;

	if (a.is_def() && b.is_def()) {
		tmp.undef_ = false;
		for (unsigned int i=0; i<7; i++) {
			tmp.exponents_[i] = a.exponents_[i] - b.exponents_[i];
		}
	}

	return tmp;
}
