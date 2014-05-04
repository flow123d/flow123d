/*
 * flag_array.hh
 *
 *  Created on: Apr 30, 2014
 *      Author: jb
 */

#ifndef FLAG_ARRAY_HH_
#define FLAG_ARRAY_HH_

#include <ostream>

/**
 * @brief std::bitset with generalized mask mechanism.
 *
 * The drawback of classical bitfield technique implemented e.g. in std::bitset<>
 * is problematic work (set and test) with more bits then one at time. This class is simple wrapper
 * for std::bitset<> that provides public class FlagArray<..>::Mask that represents some
 * bit subset together with a reference value which can be used to set the bit subset or to test the bit subset.
 *
 * The class is meant to be used in several other classes, that can define their specific masks.
 * Every class using the FlagArray should provide specific @p Tag template parameter in order to
 * guarantee Mask - FlagArray compatibility.
 *
 * Implementation note: since std::bitset do not have constexpr operators (may be in c++1y)
 * we currently use just unsigned int in our implementation.
 *
 * Usage:
 *
 * @code
 * class FlagsUser {
 * public:
 *      static constexpr
 * };
 * @endcode
 */
template <class Tag, int Size = 32>
class FlagArray {

private:
    /// Declaration of internal bitfield type.
    typedef unsigned int BitField;

public:


    /**
     * Class defines a flag mask that can set or reset certain bits to predefined values.
     * Thus it consists of a @p mask part that specifies bits which are set to the second part @p set.
     *
     * A Mask can be applied to an FlagArray via FlagArray::set() method, masks can be combined through
     * the "&" operator and can be used to test flags via FlagArray::is() method.
     */
    class Mask {
    public:

        constexpr Mask() = default;

        /// Constructor.
        constexpr Mask(BitField mask, BitField set)
        : mask_(mask), set_(set)
        {}

        /**
         * Simple constructor. Namely to get elementary mask for one bit, e.g.
         * @code
         *      static constexpr Mask input_flag{0x0010};
         * @endcode
         */
        constexpr Mask(BitField mask)
        : Mask(mask,mask)
        {}

        /**
         *  Apply mask @p other to *this mask. That is, join bit masks and overwrite
         *  bits given by @p other.mask_ by @p outher.set_ values.
         *  The action of the result Mask on the FlagArray (method @p set()) is the same as application of *this
         *  followed by the application of @p other.
         */
        constexpr Mask operator&(Mask other) const
        { return Mask( mask_ | other.mask_, mask_set(set_, other.mask_, other.set_) );
        }


        /// Mask negation.
        constexpr Mask operator~() const
        { return Mask( mask_ , mask_ & ~set_); }

        constexpr bool match(BitField flags) const
        { return ((mask_ & flags) ^ (mask_ &  set_) == 0); }

        friend std::ostream &operator<<(std::ostream &stream, const Mask &m)
        { stream << std::hex << m.mask_ << ", " << std::hex << m.set_;
          return stream;
        }

    private:

        /// Returns bitset @p set1 with bits given by @p mask set to the values given by the bitset @p set2.
        static constexpr BitField mask_set( BitField set1, BitField mask2, BitField set2)
        {  return (mask2 & set2) | (~mask2 & set1);
        }
    public:
        /// The mask operates only on true bits of the @p mask_ member.
        BitField mask_;
        /// Values of the mask. Only bits given by @p mask_ are significant.
        BitField set_;

        template <class T, int S>
        friend class FlagArray;
    };


public:
    /// Allocated size of flags storage.
    static const unsigned int size = Size;

    static constexpr Mask all_true_mask=Mask(~BitField(0), ~BitField(0));
    static constexpr Mask all_false_mask=Mask(~BitField(0), BitField(0));
    static constexpr Mask none_mask=Mask(0, 0);




    /// Default constructor turns all flags off.
    FlagArray()
    : flags_(BitField(0))
    {}

    /// Conversion from the mask.
    FlagArray(Mask mask)
    : FlagArray()
    {
        this->set(mask);
    }

    /**
     * The FlagArray match a mask if and only if
     * bits given by the @p mask.mask_ are same in both
     * the flags_ and mask.set_.
     */
    constexpr bool match(Mask mask) const
    { return mask.match(flags_); }

    /**
     * Apply the mask to the flags.
     * Bits by the @p mask.mask_ are overwritten by the @p mask.set_ values.
     */
    FlagArray &set(Mask mask)
    { flags_ = Mask::mask_set(flags_, mask.mask_, mask.set_); return *this; }

    friend std::ostream &operator<<(std::ostream & s, const FlagArray &f) {
        return (s << f.flags_);
    }

private:
    /// flags storage
    BitField flags_;
};


#endif /* FLAG_ARRAY_HH_ */
