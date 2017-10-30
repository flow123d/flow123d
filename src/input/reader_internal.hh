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
 * @file    reader_internal.hh
 * @brief
 */

#ifndef READER_INTERNAL_HH_
#define READER_INTERNAL_HH_


#include "input/input_type_forward.hh"

#include "input/storage.hh"
#include "input/path_base.hh"
#include "input/reader_internal_base.hh"


namespace Input {

using namespace std;


/**
 * @brief Creates storage of IST defined in JSON or YAML file.
 *
 * This class works like start point of creating IST storage. Allows to use other descendants
 * of ReaderInternalBase to construct storage of special parts of IST:
 *  - transposition of Type::Array type
 *  - subtree included in CSV file
 *
 * @ingroup input
 */
class ReaderInternal : public ReaderInternalBase {
public:
	/// Constructor
	ReaderInternal();

	/// Create storage of given @p type.
    StorageBase * read_storage(PathBase &p, const Type::TypeBase *type);

protected:
    StorageBase * make_sub_storage(PathBase &p, const Type::Array *array) override;           ///< Create storage of Type::Array type
    StorageBase * make_sub_storage(PathBase &p, const Type::Selection *selection) override;   ///< Create storage of Type::Selection type
    StorageBase * make_sub_storage(PathBase &p, const Type::Bool *bool_type) override;        ///< Create storage of Type::Bool type
    StorageBase * make_sub_storage(PathBase &p, const Type::Integer *int_type) override;      ///< Create storage of Type::Integer type
    StorageBase * make_sub_storage(PathBase &p, const Type::Double *double_type) override;    ///< Create storage of Type::Double type
    StorageBase * make_sub_storage(PathBase &p, const Type::String *string_type) override;    ///< Create storage of Type::String type

};


} // namespace Input



#endif /* READER_INTERNAL_HH_ */
