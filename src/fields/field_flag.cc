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
 * @file    field_flag.cc
 * @brief   
 */

#include "fields/field_flag.hh"

constexpr FieldFlag::Flags::Mask FieldFlag::equation_input;
constexpr FieldFlag::Flags::Mask FieldFlag::declare_input;
constexpr FieldFlag::Flags::Mask FieldFlag::allow_output;

constexpr FieldFlag::Flags::Mask FieldFlag::input_copy;

constexpr FieldFlag::Flags::Mask FieldFlag::in_time_term;
constexpr FieldFlag::Flags::Mask FieldFlag::in_main_matrix;
constexpr FieldFlag::Flags::Mask FieldFlag::in_rhs;

constexpr FieldFlag::Flags::Mask FieldFlag::equation_result;
constexpr FieldFlag::Flags::Mask FieldFlag::equation_external_output;
