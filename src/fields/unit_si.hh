/*
 * unit.hh
 *
 *  Created on: Aug 30, 2012
 *      Author: jb
 */


#include <vector>
#include <string>


class Unit {
public:
	/// Constructor
	Unit();

	/// Values represent symbols of base SI units in same order as units are stored in exponents_ vector
	static const std::vector<std::string> unit_symbols;

	/// Methods set values of exponents for SI units with similar name
	Unit & m(int val);
	Unit & kg(int val);
	Unit & s(int val);
	Unit & A(int val);
	Unit & K(int val);
	Unit & mol(int val);
	Unit & cd(int val);

	/**
	 * Printout units in format: m.kg^2.s^-2
	 */
	std::string print();

private:
	/**
	 * Stores exponents of base SI units in this order:
	 * [m, kg, s, A, K, mol, cd]
	 */
	std::vector<int> exponents_;

};
