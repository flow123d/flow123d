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

namespace units_converter
{

/// Store structure given by parser
struct Factor {
	/// Constructor
	Factor() : exponent_(1) {}
	Factor(std::string factor, int exponent) : factor_(factor), exponent_(exponent) {}

	std::string factor_;  //!< string represantation of unit or user defined constant
	int exponent_;        //!< exponent
};
struct Formula {
	/// Constructor
	Formula() : coef_(1.0) {}

	double coef_;                  //!< multiplicative coeficient
	std::vector<struct Factor> factors_;  //!< factors of formula
};
typedef typename std::map<std::string, struct Formula> UnitData;


}  // namespace units_converter

#endif /* UNIT_CONVERTER_HH_ */
