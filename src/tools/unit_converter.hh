/*!
 *
﻿ * Copyright (C) 2015 Technical University of Liberec.  All rights reserved.
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License version 3 as published by the
 * Free Software Foundation. (http://www.gnu.org/licenses/gpl-3.0.en.html)
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 *
 * @file    unit_converter.hh
 * @brief
 */

#ifndef UNIT_CONVERTER_HH_
#define UNIT_CONVERTER_HH_


#include <map>                       // for map, map<>::value_compare
#include <string>                    // for string
#include <vector>                    // for vector

#include "tools/unit_si.hh"         // for UnitSI
#include "input/input_exception.hh"  // for DECLARE_INPUT_EXCEPTION, Exception
#include "system/asserts.hh"         // for Assert, ASSERT
#include "system/exceptions.hh"      // for operator<<, ExcStream, EI, TYPED...

// Declaration of exceptions
TYPEDEF_ERR_INFO(EI_UnitDefinition, std::string);
TYPEDEF_ERR_INFO(EI_ExpectedUnit, std::string);
TYPEDEF_ERR_INFO(EI_UnitError, std::string);
DECLARE_INPUT_EXCEPTION(ExcInvalidUnit,
        << "Invalid definition of unit: " << EI_UnitDefinition::qval << "\n" << EI_UnitError::val << ".\n");
DECLARE_INPUT_EXCEPTION(ExcNoncorrespondingUnit,
        << "Non-corresponding definition of unit: " << EI_UnitDefinition::qval << "\nExpected: unit with base format "
		<< EI_ExpectedUnit::qval << ".\n");


/// Store structure given by parser
struct Factor {
	/// Constructor
	Factor() : exponent_(1), basic_(true) {}
	Factor(std::string factor, int exponent, bool basic = true) : factor_(factor), exponent_(exponent), basic_(basic) {}

	std::string factor_;  //!< string represantation of unit or user defined constant
	int exponent_;        //!< exponent
	bool basic_;          //!< unit is basic (strict defined in application) / derived (defined by user as formula)
};
struct Formula {
	/// Constructor
	Formula() : coef_(1.0) {}

	double coef_;                         //!< multiplicative coeficient
	std::vector<struct Factor> factors_;  //!< factors of formula
};
typedef typename std::map<std::string, struct Formula> UnitData;



/**
 * @brief Helper class. Defines basic factors of SI, non-SI and derived units.
 *
 * Class is instantiated as static object of UnitConverter class and provides check
 * of user defined units.
 */
class BasicFactors {
public:
	struct DerivedUnit {
		double coef_;    //!< multiplicative coeficient
		UnitSI unit_;    //!< derived SI unit
	};

	typedef std::map<std::string, struct DerivedUnit> UnitsMap;

	/// Define all base and derived units given by their symbol.
	UnitsMap units_map_;

	/// Constructor
	BasicFactors();

};


class UnitConverter {
public:
	/// Define all base and derived units given by their symbol.
	static const BasicFactors basic_factors;

	/// Constructor
	UnitConverter();

	/// Convert string to coeficient and UnitSI representation, return coeficient
	double convert(std::string actual_unit);

	/// Return @p unit_si_
	inline UnitSI unit_si() const {
		ASSERT(unit_si_.is_def()).error("UnitSI is not defined, first call convert method.");
		return unit_si_;
	}

protected:
	/**
	 * @brief Parse and check unit defined in string format.
	 *
	 * Return data in format \p UnitData
	 */
	UnitData read_unit(std::string s);

	/// Calculates UnitSi and coeficient of Factor, recursively calls this method for user defined formula
	void add_converted_unit(Factor factor, UnitData &unit_data, UnitSI &unit_si, double &coef);

	/**
	 * Coeficient of unit.
	 *
	 * Coeficient is used if unit is not in basic format. Example: if the unit is specified
	 * in minutes, coeficient has value 60.
	 */
	double coef_;

	/**
	 * Basic format of converted SI unit
	 */
	UnitSI unit_si_;

};


#endif /* UNIT_CONVERTER_HH_ */
