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
 * @file    reader_internal_transpose.cc
 * @brief
 */


#include "input/reader_internal_transpose.hh"
#include "input/input_type.hh"
//#include "input/accessors.hh"

namespace Input {

using namespace std;

/*******************************************************************
 * implementation of ReaderInternalTranspose
 */

ReaderInternalTranspose::ReaderInternalTranspose()
: transpose_index_(0)
{}

StorageBase * ReaderInternalTranspose::read_storage(PathBase &p, const Type::Array *array)
{
	const Type::TypeBase &sub_type = array->get_sub_type();
	StorageBase *first_item_storage;
	try {
		first_item_storage = make_storage(p, &sub_type);
	} catch (ExcInputError &e) {
		if ( !array->match_size(1) ) {
			e << EI_Specification("The value should be '" + p.get_node_type(ValueTypes::array_type) + "', but we found: ");
		}
		e << EI_TransposeIndex(transpose_index_);
		e << EI_TransposeAddress(p.as_string());
		throw;
	}

	// automatic conversion to array with one element
	if (transpose_array_sizes_.size() == 0) {
		return make_autoconversion_array_storage(p, array, first_item_storage);
	} else {

		// check that transposed arrays are of the same size
		transpose_array_sizes_.erase( unique( transpose_array_sizes_.begin(), transpose_array_sizes_.end() ),
									  transpose_array_sizes_.end() );
		if (transpose_array_sizes_.size() == 1) {
			unsigned int sizes = transpose_array_sizes_[0]; // sizes of transposed

			// array size out of bounds
			if ( !array->match_size( sizes ) ) {
				stringstream ss;
				ss << "Result of transpose auto-conversion do not fit the size " << sizes << " of the Array.";
				this->generate_input_error(p, array, ss.str(), false);
			}

			// create storage of array
			StorageArray *storage_array = new StorageArray(sizes);
			storage_array->new_item(0, first_item_storage);
			if (sizes>1) {
				++transpose_index_;
				while (transpose_index_ < sizes) {
					try {
						storage_array->new_item(transpose_index_, make_storage(p, &sub_type));
					} catch (ExcInputError &e) {
						e << EI_TransposeIndex(transpose_index_);
						e << EI_TransposeAddress(p.as_string());
						throw;
					}
					++transpose_index_;
				}
			}

			return storage_array;
		} else {
			this->generate_input_error(p, array, "Unequal sizes of sub-arrays during transpose auto-conversion of '" + p.get_node_type(ValueTypes::array_type) + "'", false);
		}
	}

	return NULL;
}

StorageBase * ReaderInternalTranspose::make_sub_storage(PathBase &p, const Type::Array *array)
{
	int arr_size;
	if ( (arr_size = p.get_array_size()) != -1 ) {
		return this->make_array_storage(p, array, arr_size);
	} else {
		// if transposition is carried, only conversion to array with one element is allowed
		// try automatic conversion to array with one element
		const Type::TypeBase &sub_type = array->get_sub_type();
		StorageBase *one_element_storage = this->make_storage(p, &sub_type);
		return make_autoconversion_array_storage(p, array, one_element_storage);
	}
}

StorageBase * ReaderInternalTranspose::make_sub_storage(PathBase &p, const Type::Selection *selection)
{
	if ( p.is_array_type() ) {
		// transpose auto-conversion for array type
		return this->make_transposed_storage(p, selection);
	} else {
	    string item_name = read_string_value(p, selection);
		try {
			int value = selection->name_to_int( item_name );
			return new StorageInt( value );
		} catch (Type::Selection::ExcSelectionKeyNotFound &exc) {
			this->generate_input_error(p, selection, "Wrong value '" + item_name + "' of the Selection.", false);
		}
	}

    return NULL;
}

StorageBase * ReaderInternalTranspose::make_sub_storage(PathBase &p, const Type::Bool *bool_type)
{
	if ( p.is_array_type() ) {
		// transpose auto-conversion for array type
		return this->make_transposed_storage(p, bool_type);
	} else {
		return new StorageBool( read_bool_value(p, bool_type) );
	}
}

StorageBase * ReaderInternalTranspose::make_sub_storage(PathBase &p, const Type::Integer *int_type)
{
	if ( p.is_array_type() ) {
		// transpose auto-conversion for array type
		return this->make_transposed_storage(p, int_type);
	}
	std::int64_t value = read_int_value(p, int_type);

	if ( int_type->match(value) )
	{
		return new StorageInt( value );
	} else {
		this->generate_input_error(p, int_type, "Value out of bounds.", false);
	}
	return NULL;
}

StorageBase * ReaderInternalTranspose::make_sub_storage(PathBase &p, const Type::Double *double_type)
{
	if ( p.is_array_type() ) {
		// transpose auto-conversion for array type
		return this->make_transposed_storage(p, double_type);
	} else {
	    double value = read_double_value(p, double_type);

	    if (double_type->match(value)) {
	        return new StorageDouble( value );
	    } else {
	    	this->generate_input_error(p, double_type, "Value out of bounds.", false);
	    }
	}

    return NULL;
}

StorageBase * ReaderInternalTranspose::make_sub_storage(PathBase &p, const Type::String *string_type)
{
	if ( p.is_array_type() ) {
		// transpose auto-conversion for array type
		return this->make_transposed_storage(p, string_type);
	} else {
		return new StorageString( read_string_value(p, string_type) );
	}
}

StorageBase * ReaderInternalTranspose::make_transposed_storage(PathBase &p, const Type::TypeBase *type) {
	ASSERT_PERMANENT(p.is_array_type()).error();

	int arr_size = p.get_array_size();
	if ( arr_size == 0 ) {
		this->generate_input_error(p, type, "Empty array during transpose auto-conversion.", false);
	} else {
		if (transpose_index_ == 0) transpose_array_sizes_.push_back( arr_size );
		p.down(transpose_index_);
		StorageBase *storage = make_storage(p, type);
		p.up();
		return storage;
	}

	return NULL;
}

StorageBase * ReaderInternalTranspose::make_autoconversion_array_storage(PathBase &p, const Type::Array *array, StorageBase *item)
{
	if ( array->match_size( 1 ) ) {
		StorageArray *storage_array = new StorageArray(1);
		storage_array->new_item(0, item);

		return storage_array;
	} else {
		this->generate_input_error(p, array, "During transpose auto-conversion, the conversion to the single element array not allowed. Require type: '" + p.get_node_type(ValueTypes::array_type) + "'\nFound on input: ", true);
	}

	return NULL;
}


} // namespace Input
