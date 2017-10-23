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

#include "input/input_exception.hh"
#include "input/storage.hh"
#include "input/path_base.hh"
#include "input/reader_to_storage.hh"


namespace Input {

using namespace std;


class ReaderInternalBase {
public:
	/// Constructor
	ReaderInternalBase();

protected:

	/**
     * @brief Create storage of given @p type.
     *
     * Check correctness of the input given by json_spirit or YAML-cpp node at head() of PathBase @p p
     * against type specification @p type. Die on input error (and return NULL).
     * For correct input, creates the storage tree and returns pointer to its root node.
     */
    StorageBase * make_storage(PathBase &p, const Type::TypeBase *type);

	virtual StorageBase * make_sub_storage(PathBase &p, const Type::Record *record);           ///< Create storage of Type::Record type
    virtual StorageBase * make_sub_storage(PathBase &p, const Type::Tuple *tuple);             ///< Create storage of Type::Tuple type
    virtual StorageBase * make_sub_storage(PathBase &p, const Type::Abstract *abstr_rec);      ///< Create storage of Type::Abstract type
    virtual StorageBase * make_sub_storage(PathBase &p, const Type::Array *array)=0;           ///< Create storage of Type::Array type
    virtual StorageBase * make_sub_storage(PathBase &p, const Type::Selection *selection)=0;   ///< Create storage of Type::Selection type
    virtual StorageBase * make_sub_storage(PathBase &p, const Type::Bool *bool_type)=0;        ///< Create storage of Type::Bool type
    virtual StorageBase * make_sub_storage(PathBase &p, const Type::Integer *int_type)=0;      ///< Create storage of Type::Integer type
    virtual StorageBase * make_sub_storage(PathBase &p, const Type::Double *double_type)=0;    ///< Create storage of Type::Double type
    virtual StorageBase * make_sub_storage(PathBase &p, const Type::String *string_type)=0;    ///< Create storage of Type::String type

    /// Apply automatic conversion of Type::Record type
    StorageBase * record_automatic_conversion(PathBase &p, const Type::Record *record);

    /// Apply automatic conversion of Type::Abstract type
    StorageBase * abstract_automatic_conversion(PathBase &p, const Type::Abstract *abstr_rec);

    /// Create storage of Type::Array with given size
    StorageBase * make_array_storage(PathBase &p, const Type::Array *array, int arr_size);

    /// Dispatch according to @p type and create corresponding storage from the given string.
    StorageBase * make_storage_from_default( const string &dflt_str, std::shared_ptr<Type::TypeBase> type);

    /// Create storage of included YAML or JSON input file
    StorageBase * make_include_storage(PathBase &p, const Type::Record *record);

    bool read_bool_value(PathBase &p, const Type::TypeBase *type);           ///< Read boolean value from path
    std::int64_t read_int_value(PathBase &p, const Type::TypeBase *type);    ///< Read integer value from path
    double read_double_value(PathBase &p, const Type::TypeBase *type);       ///< Read double value from path
    std::string read_string_value(PathBase &p, const Type::TypeBase *type);  ///< Read string value from path

    /// Helper method. Get string value of included file or throw exception if reading failed.
    std::string get_included_file(PathBase &p);

};


class ReaderInternal : public ReaderInternalBase {
public:
	/// Constructor
	ReaderInternal();

protected:
    StorageBase * make_sub_storage(PathBase &p, const Type::Array *array) override;           ///< Create storage of Type::Array type
    StorageBase * make_sub_storage(PathBase &p, const Type::Selection *selection) override;   ///< Create storage of Type::Selection type
    StorageBase * make_sub_storage(PathBase &p, const Type::Bool *bool_type) override;        ///< Create storage of Type::Bool type
    StorageBase * make_sub_storage(PathBase &p, const Type::Integer *int_type) override;      ///< Create storage of Type::Integer type
    StorageBase * make_sub_storage(PathBase &p, const Type::Double *double_type) override;    ///< Create storage of Type::Double type
    StorageBase * make_sub_storage(PathBase &p, const Type::String *string_type) override;    ///< Create storage of Type::String type

};


class ReaderInternalTranspose : public ReaderInternalBase {
public:
	/// Constructor
	ReaderInternalTranspose();

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

};


class ReaderInternalCsvInclude : public ReaderInternalBase {
public:
	/// Constructor
	ReaderInternalCsvInclude();

protected:
    StorageBase * make_sub_storage(PathBase &p, const Type::Array *array) override;           ///< Create storage of Type::Array type
    StorageBase * make_sub_storage(PathBase &p, const Type::Selection *selection) override;   ///< Create storage of Type::Selection type
    StorageBase * make_sub_storage(PathBase &p, const Type::Bool *bool_type) override;        ///< Create storage of Type::Bool type
    StorageBase * make_sub_storage(PathBase &p, const Type::Integer *int_type) override;      ///< Create storage of Type::Integer type
    StorageBase * make_sub_storage(PathBase &p, const Type::Double *double_type) override;    ///< Create storage of Type::Double type
    StorageBase * make_sub_storage(PathBase &p, const Type::String *string_type) override;    ///< Create storage of Type::String type

    /// Create vector which contains actual indexes of subtree imported in CSV file.
    vector<unsigned int> create_indexes_vector(PathBase &p);

    /// Depth of CSV included subtree
    unsigned int csv_subtree_depth_;

    /// Helper vector which contains actual indexes of subtree imported in CSV file.
    vector<unsigned int> csv_storage_indexes_;

    /// Map of columns in CSV file to storage of subtree
    map<unsigned int, ReaderToStorage::IncludeCsvData> csv_columns_map_;

};


} // namespace Input



#endif /* READER_INTERNAL_HH_ */
