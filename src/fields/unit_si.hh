/*
 * unit_si.hh
 *
 *  Created on: Aug 30, 2012
 *      Author: jb
 */


#ifndef UNIT_SI_HH_
#define UNIT_SI_HH_

#include <vector>
#include <string>


class UnitSI {
public:
	/// Constructor
	UnitSI();

	/// Definitions of frequently used derived units
	/// Definition of Newton
	static const UnitSI N;
	/// Definition of Joule
	static const UnitSI J;
	/// Definition of Watt
	static const UnitSI W;
	/// Definition of Pascal
	static const UnitSI Pa;
	/// Definition of dimensionless unit
	static const UnitSI dimensionless;

	/// Values represent symbols of base SI units in same order as units are stored in exponents_ vector
	static const std::vector<std::string> unit_symbols;

	/// Methods set values of exponents for SI units with similar name
	UnitSI & m(int exp = 1);
	UnitSI & kg(int exp = 1);
	UnitSI & s(int exp = 1);
	UnitSI & A(int exp = 1);
	UnitSI & K(int exp = 1);
	UnitSI & mol(int exp = 1);
	UnitSI & cd(int exp = 1);

	/**
	 * Makes unit description string in Latex format, e.g. "$[m.kg^{2}.s^{-2}]$"
	 *
	 * Have assert for undefined units.
	 */
	std::string format() const;

	/// Set flag that unit is undefined.
	void undef(bool val = true);

	/// Return true if the unit is defined.
	bool is_def() const;

private:
	/**
	 * Stores exponents of base SI units in this order:
	 * [m, kg, s, A, K, mol, cd]
	 */
	std::vector<int> exponents_;

	/**
	 * Flag if object is undefined.
	 *
	 * Value is set on true in constructor, when any exponent is changed, false value is set.
	 */
	bool undef_;

	/// Product of two units.
	friend UnitSI operator *(const UnitSI &a, const UnitSI &b);
	/// Proportion of two units.
	friend UnitSI operator /(const UnitSI &a, const UnitSI &b);
};


/// Product of two units.
UnitSI operator *(const UnitSI &a, const UnitSI &b);

/// Proportion of two units.
UnitSI operator /(const UnitSI &a, const UnitSI &b);


#endif /* UNIT_SI_HH_ */
