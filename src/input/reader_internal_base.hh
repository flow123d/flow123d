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
 * @file    reader_internal_base.hh
 * @brief
 */

#ifndef READER_INTERNAL_BASE_HH_
#define READER_INTERNAL_BASE_HH_


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
    DECLARE_INPUT_EXCEPTION( ExcDuplicitTag, << "Error in input file: " << EI_File::qval << " at address: " << EI_ErrorAddress::qval << "\n"
                                              << "Duplicit Tag: " << EI_Tag::qval << "\n" );
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
    StorageBase * make_include_storage(PathBase &p, const Type::TypeBase *type);

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


} // namespace Input

#endif /* READER_INTERNAL_BASE_HH_ */
