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
 * @file    field_flag.hh
 * @brief   
 */

#ifndef FIELD_FLAGS_HH_
#define FIELD_FLAGS_HH_

#include "system/flag_array.hh"

class FieldFlag
{
public:
    typedef FlagArray<FieldFlag> Flags;
    typedef Flags::Mask Mask;

    /// Number of bits used by Field itself.
    static constexpr unsigned int flags_size_ = 3;

    /// The field is data parameter of the owning equation. (default on)
    static constexpr Mask equation_input{1 << 0};
    /// The field can be set from input. The key in input field descriptor is declared. (default on)
    static constexpr Mask declare_input{1 << 1};
    /// The field can output. Is part of generated output selection. (default on)
    static constexpr Mask allow_output{1 << 2};


    /**
     *  A field that is input of its equation and can not read from input, thus must be set by copy.
     *  Disabled flag 'declare_input' still allow initialization of the field from a default value.
     */
    static constexpr Mask input_copy = ~declare_input & equation_input;

    /// A field is part of time term of the equation.
    static constexpr Mask in_time_term{1 << 8};
    /// A field is part of main "stiffness matrix" of the equation.
    static constexpr Mask in_main_matrix{1 << 9};
    /// A field is part of the right hand side of the equation.
    static constexpr Mask in_rhs{1 << 10};


    /// Match result fields. These are never given by input or copy of input.
    static constexpr Mask equation_result =  allow_output & ~declare_input & ~equation_input;

    /// Match an output field, that can be also copy of other field.
    static constexpr Mask equation_external_output =  allow_output & input_copy;

};







#endif /* FIELD_FLAGS_HH_ */
