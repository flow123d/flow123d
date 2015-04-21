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


/**
 * @brief Class for representation SI units of Fields.
 *
 * Units are set through exponents of basic SI units. These exponents are set in methods with same
 * name as unit symbols (e.g. kg(), K() etc).
 *
 * Class contains method that provides formated string representing full unit symbol (usable in
 * LaTeX output).
 *
 * UnitSI object contains flag that says if it is defined
 * - undefined object can't be formated
 * - if any exponent is set, flag is set to defined
 *
 * Class contains static methods that return frequently used derived units (Watt, Pascal etc).
 */
class UnitSI {
public:
	/// Constructor
	UnitSI();

	/// Methods return frequently used derived units
	/// Returns Newton
	static UnitSI & N();
	/// Returns Joule
	static UnitSI & J();
	/// Returns Watt
	static UnitSI & W();
	/// Returns Pascal
	static UnitSI & Pa();
	/// Returns dimensionless unit
	static UnitSI & dimensionless();

	/// Methods set values of exponents for SI units with similar name
	UnitSI & m(int exp = 1);
	UnitSI & kg(int exp = 1);
	UnitSI & s(int exp = 1);
	UnitSI & A(int exp = 1);
	UnitSI & K(int exp = 1);
	UnitSI & mol(int exp = 1);
	UnitSI & cd(int exp = 1);
	/// Method sets value of exponent for m^{-d}, where d is dimension of region
	UnitSI & md(int exp = -1);

	/**
	 * Makes unit description string in Latex format, e.g. "$[m.kg^{2}.s^{-2}]$"
	 *
	 * Have assert for undefined units.
	 */
	std::string format_latex() const;

	std::string format_text() const;

	/**
	 * Set flag that unit is undefined.
	 *
	 * Default value is true (set in constructor).
	 * If any exponent is set, @p undef_ flag is unset.
	 * In all fields unit must be defined by user.
	 */
	void undef(bool val = true);

	/// Return true if the unit is defined.
	bool is_def() const;

private:
	/// Values determine positions of exponents in exponents_ vector
	enum UnitOrder {
		order_m=0,
		order_md=1,
		order_kg=2,
		order_s=3,
		order_A=4,
		order_K=5,
		order_mol=6,
		order_cd=7,
		n_base_units=8
	};

	/// Variable parts of output format. Used in the @p format method.
	struct OutputFormat {
	    std::string exp_open, exp_close, delimiter;
	};

	/// Generic output formating mtehod.
	std::string format(OutputFormat form) const;

	/**
	 * Stores exponents of base SI units in this order:
	 * [m, kg, s, A, K, mol, cd, md]
	 *
	 * where md represents value of exponent depended on dimension (m^{-d})
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
