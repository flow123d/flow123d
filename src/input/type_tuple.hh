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
 * @file    type_tuple.hh
 * @brief
 */

#ifndef TYPE_TUPLE_HH_
#define TYPE_TUPLE_HH_

#include "type_record.hh"

namespace Input {

namespace Type {

/**
 * @brief Tuple type proxy class.
 *
 * Class is used for small data items, which are repeated several times. Declaration is same
 * as for the Record. Allows two possibilities of declaration:
 *  - same as the Record notation (using named keys)
 *  - list of values in the order how the keys are defined
 *
 * In declaration, obligatory keys must be given as first, optional and default at read time
 * keys follows. Correct order is checked in @p finish method. Auto conversion is allowed
 * for the first key.
 *
 * Tuple type can't be descendant of Abstract.
 */
class Tuple : public Record {
public:
    /*
     * Exceptions specific to this class.
     */
    TYPEDEF_ERR_INFO( EI_TupleName, const string);
    DECLARE_EXCEPTION( ExcTupleWrongKeysOrder, << "Incorrect order of obligatory and non-obligatory keys in Tuple: " <<  EI_TupleName::qval );


    /// Default constructor. Empty handle.
	Tuple();

    /**
     * @brief Copy constructor.
     *
     * We allow only copies of non-empty tuples.
     */
	Tuple(const Tuple & other);


    /**
     * @brief Basic constructor.
     *
     * You have to provide \p type_name of the new declared Record type and its \p description.
     */
	Tuple(const string & type_name_in, const string & description);

    /**
     * Implements @p TypeBase::content_hash.
     *
     * Hash is calculated by type name, description, hash of keys and attributes.
     */
    TypeHash content_hash() const  override;

	/// Override @p Type::TypeBase::class_name.
	string class_name() const override { return "Tuple"; }

	/**
	 * @brief Override Record::allow_auto_conversion
	 *
	 * Forbids setting of auto conversion key. Tuple type can have auto convertible the first key.
	 */
	Tuple &allow_auto_conversion(const string &from_key) override;

    /// Class comparison and Tuple type name comparison.
    bool operator==(const TypeBase &other) const;

    /**
     * @brief Close the Tuple for further declarations of keys.
     *
     * @see Record::close
     */
    const Tuple &close() const;

    /**
     * @brief Finish declaration of the Tuple type.
     *
     * Completes Tuple (check auto convertible key, parameters of generic types etc).
     */
    bool finish(bool is_generic = false) override;

    /**
     * @brief Override Record::derive_from
     *
     * Deriving of Tuple type is forbidden. Type is determined for small simple data.
     */
    Tuple &derive_from(Abstract &parent);

    /**
     * @brief Return count of obligatory keys.
     *
     * Needs in exceptions.
     */
    unsigned int obligatory_keys_count() const;

    /// Implements @p TypeBase::make_instance.
    MakeInstanceReturnType make_instance(std::vector<ParameterPair> vec = std::vector<ParameterPair>()) const override;

    /// Overrides Record::deep_copy
    Tuple deep_copy() const;

    /// Overrides Record::root_of_generic_subtree
    Tuple &root_of_generic_subtree() override;

    /// Overrides Record::declare_key(const string &, boost::shared_ptr<TypeBase>, const Default &, const string &)
    Tuple &declare_key(const string &key, boost::shared_ptr<TypeBase> type,
                            const Default &default_value, const string &description);

    /// Overrides Record::declare_key(const string &, const KeyType &, const Default &, const string &)
    template <class KeyType>
    Tuple &declare_key(const string &key, const KeyType &type,
                            const Default &default_value, const string &description);


    /// Overrides Record::declare_key(const string &, const KeyType &, const string &)
    template <class KeyType>
    Tuple &declare_key(const string &key, const KeyType &type,
                            const string &description);

};

} // closing namespace Type
} // closing namespace Input

#endif /* TYPE_TUPLE_HH_ */
