/*
 * field_flags.hh
 *
 *  Created on: May 3, 2014
 *      Author: jb
 */

#ifndef FIELD_FLAGS_HH_
#define FIELD_FLAGS_HH_

#include "system/flag_array.hh"

class FieldFlag
{
public:
    typedef FlagArray<FieldFlag> Flags;

    /// Number of bits used by Field itself.
    static constexpr unsigned int flags_size_ = 3;

    /// The field is data parameter of the owning equation. (default on)
    static constexpr Flags::Mask equation_input{1};
    /// The field can be set from input. The key in input field descriptor is declared. (default on)
    static constexpr Flags::Mask declare_input{2};
    /// The field can output. Is part of generated output selection. (default on)
    static constexpr Flags::Mask allow_output{4};

    static constexpr Flags::Mask in_time_term{8};
    static constexpr Flags::Mask in_mass_matrix{16};
    static constexpr Flags::Mask in_right_hand_side{32};

    /// Match non-result fields, that are data fields of an equation.
    static constexpr Flags::Mask equation_result =  ~declare_input & ~equation_input;
};







#endif /* FIELD_FLAGS_HH_ */
