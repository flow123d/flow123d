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
 * @file    unit_converter.cc
 * @brief
 */


#include "fields/unit_converter.hh"

const std::map<std::string, UnitConverter::DerivedUnit> UnitConverter::units_map = {
		{ "mm",  { 0.001,   UnitSI().m() } },
		{ "cm",  { 0.01,    UnitSI().m() } },
		{ "dm",  { 0.1,     UnitSI().m() } },
		{ "m",   { 1,       UnitSI().m() } },
		{ "km",  { 1000,    UnitSI().m() } },

		{ "g",   { 0.001,   UnitSI().kg() } },
		{ "kg",  { 1,       UnitSI().kg() } },
		{ "t",   { 1000,    UnitSI().kg() } },

		{ "s",   { 1,       UnitSI().s() } },
		{ "min", { 60,      UnitSI().s() } },
		{ "h",   { 3600,    UnitSI().s() } },
		{ "d",   { 24*3600, UnitSI().s() } },

		{ "A",   { 1,       UnitSI().A() } },

		{ "K",   { 1,       UnitSI().K() } },

		{ "cd",  { 1,       UnitSI().cd() } },

		{ "mol", { 1,       UnitSI().mol() } },

		{ "N",   { 1,       UnitSI().m().kg().s(-2) } },
		{ "kN",  { 1000,    UnitSI().m().kg().s(-2) } },

		{ "J",   { 1,       UnitSI().m(2).kg().s(-2) } },
		{ "kJ",  { 1000,    UnitSI().m(2).kg().s(-2) } },
		{ "MJ",  { 1000000, UnitSI().m(2).kg().s(-2) } },

		{ "W",   { 1,       UnitSI().m(2).kg().s(-3) } },
		{ "kW",  { 1000,    UnitSI().m(2).kg().s(-3) } },
		{ "MW",  { 1000000, UnitSI().m(2).kg().s(-3) } },

		{ "Pa",  { 1,       UnitSI().m(-1).kg().s(-2) } },
		{ "hPa", { 100,     UnitSI().m(-1).kg().s(-2) } },
		{ "kPa", { 1000,    UnitSI().m(-1).kg().s(-2) } },
		{ "MPa", { 1000000, UnitSI().m(-1).kg().s(-2) } },

		{ "rad", { 1,       UnitSI().m(0) } }
};


UnitConverter::UnitConverter(UnitSI unit_si)
: coef_(1.0), unit_si_(unit_si) {}

