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


/**
 *  @brief Internal part of IST reader.
 *
 *  This class defines interface for creating storage. Dispatching of storage is ensured by the non-virtual
 *  method read_storage, different for individual descendants.
 *
 *  @ingroup input
 */
class ReaderInternalBase {
public:
	/*
	 * Exceptions.
	 */
	TYPEDEF_ERR_INFO(EI_InputType, string );
	TYPEDEF_ERR_INFO(EI_File, const string);
	TYPEDEF_ERR_INFO(EI_Specification, const string);
	TYPEDEF_ERR_INFO(EI_Format, const string);
	TYPEDEF_ERR_INFO(EI_JSON_Type, const string);
	TYPEDEF_ERR_INFO(EI_ErrorAddress, string);
    TYPEDEF_ERR_INFO(EI_TransposeIndex, unsigned int);
    TYPEDEF_ERR_INFO(EI_TransposeAddress, string);
    TYPEDEF_ERR_INFO(EI_JSONLine, unsigned int);
    TYPEDEF_ERR_INFO(EI_JSONColumn, unsigned int);
    TYPEDEF_ERR_INFO(EI_JSONReason, string);
    TYPEDEF_ERR_INFO(EI_InputErrorMessage, const string);
    TYPEDEF_ERR_INFO(EI_RecordName, const string);
    TYPEDEF_ERR_INFO(EI_Tag, string);
	TYPEDEF_ERR_INFO(EI_TokenizerMsg, std::string);
	TYPEDEF_ERR_INFO(EI_ColumnIndex, unsigned int);

	/// General exception during conversion from JSON/YAML tree to storage.
	DECLARE_INPUT_EXCEPTION( ExcInputError, << "Error in input file: " << EI_File::qval << " at address: " << EI_ErrorAddress::qval << "\n"
	                                        << EI_Specification::val << "\n"
	                                        << EI_Format::val << " type: " << EI_JSON_Type::qval << "\n"
	                                        << "Expected type:\n" << EI_InputType::val );

	DECLARE_INPUT_EXCEPTION( ExcNotJSONFormat, << "Not valid JSON file " << EI_File::qval << ". Error at line "
            << EI_JSONLine::val << " : col " << EI_JSONColumn::val
            << " ; reason: " << EI_JSONReason::val << "\n" );
    DECLARE_INPUT_EXCEPTION( ExcAutomaticConversionError, << "Error during automatic conversion of "
    		<< EI_RecordName::val << " record.\n " << EI_InputErrorMessage::val << "\n" );
    DECLARE_INPUT_EXCEPTION( ExcForbiddenTag, << "Tag " << EI_Tag::qval << " " << EI_Specification::val << "\n" );
    DECLARE_INPUT_EXCEPTION( ExcWrongCsvFormat, << EI_Specification::val << ",\n" << EI_TokenizerMsg::val << "\n" );
    DECLARE_INPUT_EXCEPTION( ExcMultipleDefinitionCsvColumn, << "Multiple definition of column with index " << EI_ColumnIndex::qval
    		<< " in included CSV file:\n" << EI_File::val << ",\n" );

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

    /// Create storage of Type::Record type. Common method of all descendants.
	StorageBase * make_sub_storage(PathBase &p, const Type::Record *record);

	/// Create storage of Type::Tuple type. Common method of all descendants.
    StorageBase * make_sub_storage(PathBase &p, const Type::Tuple *tuple);

    /// Create storage of Type::Abstract type. Common method of all descendants.
    StorageBase * make_sub_storage(PathBase &p, const Type::Abstract *abstr_rec);

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

    /// Generate @p ExcInputError
    void generate_input_error(PathBase &p, const Type::TypeBase *type, std::string spec, bool add_type);

    /// Complete specification, error address and JSON type error tags to @p ExcInputError
    void complete_input_error(ExcInputError & e, PathBase &p, ValueTypes value_type);

};


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


/**
 * @brief Creates storage of part of subtree defined in CSV file.
 *
 * We use this class if input tree contains included CSV file where Array
 * is expected. This type must correspond with type_of_value of Array of Records.
 *
 * @ingroup input
 */
class ReaderInternalCsvInclude : public ReaderInternalBase {
public:
	/// Constructor
	ReaderInternalCsvInclude();

	/**
	 * Create storage of subtree defined in CSV file of given @p array.
	 *
	 * Processing of subtree included CSV file:
	 * 1. We expect Record with special structure and tag 'include_csv' in place of Array.
	 * 2. We load path to CSV file and number of head lines to skip in this file.
	 * 3. We load Record 'format' which determines positions of columns in CSV file.
	 *    Structure of this Record must be equal with subtype of array.
	 * 4. We fill helper @p csv_columns_map_ which allows mapping between columns in CSV
	 *    file and storage of array subtype. We return helper storage of this subtype with
	 *    default values (only determines structure).
	 * 5. We iterate through CSV file. For every line we fill data into helper storage and
	 *    copy this helper storage to storage represents included array.
	 */
    StorageBase * read_storage(PathBase &p, const Type::Array *array);

protected:
    /// List of data types, used for mapping columns in CSV include.
    typedef enum {
    	type_int, type_double, type_bool, type_string, type_sel
    } IncludeDataTypes;

    /// Data of one column of including CSV file.
    struct IncludeCsvData {
    	IncludeDataTypes data_type;
    	vector<unsigned int> storage_indexes;
    	const Type::TypeBase *type;
    };

    StorageBase * make_sub_storage(PathBase &p, const Type::Array *array) override;           ///< Create storage of Type::Array type
    StorageBase * make_sub_storage(PathBase &p, const Type::Selection *selection) override;   ///< Create storage of Type::Selection type
    StorageBase * make_sub_storage(PathBase &p, const Type::Bool *bool_type) override;        ///< Create storage of Type::Bool type
    StorageBase * make_sub_storage(PathBase &p, const Type::Integer *int_type) override;      ///< Create storage of Type::Integer type
    StorageBase * make_sub_storage(PathBase &p, const Type::Double *double_type) override;    ///< Create storage of Type::Double type
    StorageBase * make_sub_storage(PathBase &p, const Type::String *string_type) override;    ///< Create storage of Type::String type

    /// Create vector which contains actual indexes of subtree imported in CSV file.
    vector<unsigned int> create_indexes_vector(PathBase &p);

    /// Set storage of simple input type with value given from CSV file.
    void set_storage_from_csv(unsigned int column_index, StorageBase * item_storage, StorageBase * new_storage);

    /// Depth of CSV included subtree
    unsigned int csv_subtree_depth_;

    /// Map of columns in CSV file to storage of subtree
    map<unsigned int, IncludeCsvData> csv_columns_map_;

};


} // namespace Input



#endif /* READER_INTERNAL_HH_ */
