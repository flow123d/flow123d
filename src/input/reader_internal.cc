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
 * @file    reader_internal.cc
 * @brief
 */


#include "input/reader_internal.hh"
#include "input/reader_internal_transpose.hh"
#include "input/reader_internal_csv.hh"
#include "input/input_type.hh"

namespace Input {

using namespace std;


/*******************************************************************
 * implementation of ReaderInternal
 */

ReaderInternal::ReaderInternal()
{}

StorageBase * ReaderInternal::read_storage(PathBase &p, const Type::TypeBase *type)
{
	return this->make_storage(p, type);
}

StorageBase * ReaderInternal::make_sub_storage(PathBase &p, const Type::Array *array)
{
	int arr_size;
	if ( (arr_size = p.get_array_size()) != -1 ) {
		return this->make_array_storage(p, array, arr_size);
	} else if (p.get_record_tag() == "include") {
		return make_include_storage(p, array);
	} else if (p.get_record_tag() == "include_csv") {
		ReaderInternalCsvInclude reader_csv;
		return reader_csv.read_storage(p, array);
	} else {
		ReaderInternalTranspose reader_transpose;
		return reader_transpose.read_storage(p, array);
	}
}

StorageBase * ReaderInternal::make_sub_storage(PathBase &p, const Type::Selection *selection)
{
    string item_name = read_string_value(p, selection);
	try {
		int value = selection->name_to_int( item_name );
		return new StorageInt( value );
	} catch (Type::Selection::ExcSelectionKeyNotFound &exc) {
		this->generate_input_error(p, selection, "Wrong value '" + item_name + "' of the Selection.", false);
	}

    return NULL;
}

StorageBase * ReaderInternal::make_sub_storage(PathBase &p, const Type::Bool *bool_type)
{
	return new StorageBool( read_bool_value(p, bool_type) );
}

StorageBase * ReaderInternal::make_sub_storage(PathBase &p, const Type::Integer *int_type)
{
	std::int64_t value = read_int_value(p, int_type);

	if ( int_type->match(value) )
	{
		return new StorageInt( value );
	} else {
		this->generate_input_error(p, int_type, "Value out of bounds.", false);
	}
	return NULL;
}

StorageBase * ReaderInternal::make_sub_storage(PathBase &p, const Type::Double *double_type)
{
    double value = read_double_value(p, double_type);

    if (double_type->match(value)) {
        return new StorageDouble( value );
    } else {
    	this->generate_input_error(p, double_type, "Value out of bounds.", false);
    }

    return NULL;
}

StorageBase * ReaderInternal::make_sub_storage(PathBase &p, const Type::String *string_type)
{
	string value = read_string_value(p, string_type);

	if (string_type->match(value))
		return new StorageString( value );
	else
		THROW( ExcInputError() << EI_Specification("Output file can not be given by absolute path: '" + value + "'")
						<< EI_ErrorAddress(p.as_string()) << EI_JSON_Type("") << EI_InputType(string_type->desc()) );

	return NULL;
}


} // namespace Input
