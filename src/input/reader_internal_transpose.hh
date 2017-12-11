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
 * @file    reader_internal_transpose.hh
 * @brief
 */

#ifndef READER_INTERNAL_TRANSPOSE_HH_
#define READER_INTERNAL_TRANSPOSE_HH_


#include "input/input_type_forward.hh"

#include "input/storage.hh"
#include "input/path_base.hh"
#include "input/reader_internal_base.hh"


namespace Input {

using namespace std;

/**
 * @brief Creates storage of transposed subtree defined on input.
 *
 * We use this class if input tree contains another type at position where Array
 * is expected. This type must correspond with type_of_value of Array.
 *
 * @ingroup input
 */
class ReaderInternalTranspose : public ReaderInternalBase {
public:
	/// Constructor
	ReaderInternalTranspose();

	/**
	 * @brief Create storage of transposed subtree of given @p Array.
	 *
	 * Processing of subtree with transposition:
	 * 1. We set @p transpose_index_ to value '0' (transposition of first Array item).
	 * 2. We retrieve whole subtree and find Array types that are located at position
	 *    where other type is expected (type_of_value of found Array must corresponds
	 *    with excepted type).
	 *    We create storage corresponding with subtree (unexpected Arrays are replaced
	 *    by item at position given by @p transpose_index_.
	 * 3. Together with paragraph 2 we store sizes of found Arrays to
	 *    @p transpose_array_sizes_.
	 * 4. We check sizes stored in transpose_array_sizes_ (all must be in equal
	 *    and may not be equal to zero). This size determines size of transposed Array
	 *    type.
	 * 5. We repeat paragraph 2 for all items of transposed Array (gradual increase of
	 *    @p transpose_index_).
	 */
    StorageBase * read_storage(PathBase &p, const Type::Array *array);

protected:
    StorageBase * make_sub_storage(PathBase &p, const Type::Array *array) override;           ///< Create storage of Type::Array type
    StorageBase * make_sub_storage(PathBase &p, const Type::Selection *selection) override;   ///< Create storage of Type::Selection type
    StorageBase * make_sub_storage(PathBase &p, const Type::Bool *bool_type) override;        ///< Create storage of Type::Bool type
    StorageBase * make_sub_storage(PathBase &p, const Type::Integer *int_type) override;      ///< Create storage of Type::Integer type
    StorageBase * make_sub_storage(PathBase &p, const Type::Double *double_type) override;    ///< Create storage of Type::Double type
    StorageBase * make_sub_storage(PathBase &p, const Type::String *string_type) override;    ///< Create storage of Type::String type

    /// Apply transposition and create storage of Type::Array type
    StorageBase * make_transposed_storage(PathBase &p, const Type::TypeBase *type);

    /// Apply conversion to one element storage of Type::Array type
    StorageBase * make_autoconversion_array_storage(PathBase &p, const Type::Array *array, StorageBase *item);

    /// Index of processed item in transposed part of input tree.
    unsigned int transpose_index_;

    /// Helper vector what allows check sizes of all transposed Arrays.
    vector<unsigned int> transpose_array_sizes_;

};



} // namespace Input

#endif /* READER_INTERNAL_TRANSPOSE_HH_ */
