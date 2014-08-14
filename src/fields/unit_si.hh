/*
 * unit_si.hh
 *
 *  Created on: Aug 30, 2012
 *      Author: jb
 */


#include <vector>
#include <string>


class UnitSI {
public:
	/// Constructor
	UnitSI();

	/// Definitions of frequently used derived units (Newton, Joule, Watt, Pascal)
	static const UnitSI N;
	static const UnitSI J;
	static const UnitSI W;
	static const UnitSI Pa;

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
	 * Printout units in format: $[m.kg^{2}.s^{-2}]$
	 */
	std::string print() const;

	/// Set flag that object is undefined.
	void undef(bool val = true);

	// Return value of defined_ parameter.
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

friend UnitSI operator *(const UnitSI &a, const UnitSI &b);
};


UnitSI operator *(const UnitSI &a, const UnitSI &b);
