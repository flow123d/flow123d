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

#include <vector>
#include <string>
#include <map>

#include "input/input_exception.hh"
#include "fields/unit_si.hh"

namespace units_converter
{
// Declaration of exception
TYPEDEF_ERR_INFO(EI_UnitDefinition, std::string);
TYPEDEF_ERR_INFO(EI_UnitError, std::string);
DECLARE_INPUT_EXCEPTION(ExcInvalidUnit,
        << "Invalid definition of unit: " << EI_UnitDefinition::qval << "\n" << EI_UnitError::val << ".\n");


/// Store structure given by parser
struct Factor {
	/// Constructor
	Factor() : exponent_(1), basic_(true) {}
	Factor(std::string factor, int exponent) : factor_(factor), exponent_(exponent), basic_(true) {}

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

}  // namespace units_converter


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

	/// Define all base and derived units given by their symbol.
	std::map<std::string, DerivedUnit> units_map_;

	// Constructor
	BasicFactors();

};


class UnitConverter {
public:
	struct DerivedUnit {
		double coef_;    //!< multiplicative coeficient
		UnitSI unit_;    //!< derived SI unit
	};

	/// Define all base and derived units given by their symbol.
	static const BasicFactors basic_factors;

	// Constructor
	UnitConverter(UnitSI unit_si);

private:
	/**
	 * Coeficient of unit.
	 *
	 * Coeficient is used if unit is not in basic format. Example: if the unit is specified
	 * in minutes, coeficient has value 60.
	 */
	double coef_;

	/**
	 * Basic format of SI unit
	 */
	UnitSI unit_si_;

};


#endif /* UNIT_CONVERTER_HH_ */
