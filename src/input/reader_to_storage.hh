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
 * @file    reader_to_storage.hh
 * @brief   
 */

#ifndef READER_TO_STORAGE_HH_
#define READER_TO_STORAGE_HH_


#include <sstream>

#include "input/input_type_forward.hh"

#include "input/input_exception.hh"
#include "input/storage.hh"
#include "input/path_base.hh"

#include "system/file_path.hh"


namespace Input {

using namespace std;

class ReaderInternalBase;
class ReaderInternalCsvInclude;



/// Possible formats of input files.
typedef enum  {
    format_JSON,
    format_YAML
} FileFormat;



/**
 * @brief Enum of possible input types.
 *
 * Values in @p json_type_names must be stored in same order.
 */
typedef enum {
	obj_type, array_type, str_type, bool_type, int_type, real_type, null_type, scalar_type, undef_type
} ValueTypes;



/**
 *  @brief Reader for (slightly) modified input files.
 *
 *  This class implements a reader of modified input file format (JSON or YAML). The modifications include
 *  shell-like comments (using hash '#' character), this is implemented in comment_filter.hh,  optional quoting of
 *  keys in JSON objects that do not contain spaces, and possibility to use '=' instead of ':'. So you can write:
 *  @code
 *    { key1="text", key2=2, "key 3"=4 }
 *  @endcode
 *  Note, however, that our input interface allows only C identifiers for keys. The reader use json_spirit library
 *  (based on Spirit parser from Boost) with slightly modified grammar.
 *
 *  The input file is at first read and parsed by json_spirit. Then ReaderToStorage pass through tree with parsed data along
 *  with passing through declaration tree. The input data are check against declaration and stored in the Storage tree.
 *
 *  Accessor to the root record is provided by ReaderToStorage::get_root_interface<T> method template.
 *
 *  @ingroup input
 */
class ReaderToStorage {
public:
    /*
     * Exceptions.
     */
    /// General exception during conversion from JSON tree to storage.
    TYPEDEF_ERR_INFO(EI_InputType, string );
    TYPEDEF_ERR_INFO(EI_File, const string);
    TYPEDEF_ERR_INFO(EI_Specification, const string);
    TYPEDEF_ERR_INFO(EI_Format, const string);
    TYPEDEF_ERR_INFO(EI_JSON_Type, const string);
    TYPEDEF_ERR_INFO( EI_ErrorAddress, string);
    TYPEDEF_ERR_INFO( EI_TransposeIndex, unsigned int);
    TYPEDEF_ERR_INFO( EI_TransposeAddress, string);
    DECLARE_INPUT_EXCEPTION( ExcInputError, << "Error in input file: " << EI_File::qval << " at address: " << EI_ErrorAddress::qval << "\n"
                                            << EI_Specification::val << "\n"
                                            << EI_Format::val << " type: " << EI_JSON_Type::qval << "\n"
                                            << "Expected type:\n" << EI_InputType::val );


    TYPEDEF_ERR_INFO( EI_JSONLine, unsigned int);
    TYPEDEF_ERR_INFO( EI_JSONColumn, unsigned int);
    TYPEDEF_ERR_INFO( EI_JSONReason, string);
    DECLARE_INPUT_EXCEPTION( ExcNotJSONFormat, << "Not valid JSON file " << EI_File::qval << ". Error at line "
            << EI_JSONLine::val << " : col " << EI_JSONColumn::val
            << " ; reason: " << EI_JSONReason::val << "\n" );

    TYPEDEF_ERR_INFO( EI_InputErrorMessage, const string);
    TYPEDEF_ERR_INFO( EI_RecordName, const string);
    DECLARE_INPUT_EXCEPTION( ExcAutomaticConversionError, << "Error during automatic conversion of "
    		<< EI_RecordName::val << " record.\n " << EI_InputErrorMessage::val << "\n" );

    TYPEDEF_ERR_INFO( EI_Tag, string);
    DECLARE_INPUT_EXCEPTION( ExcForbiddenTag, << "Tag " << EI_Tag::qval << " " << EI_Specification::val << "\n" );

	TYPEDEF_ERR_INFO(EI_TokenizerMsg, std::string);
    DECLARE_INPUT_EXCEPTION(ExcWrongCsvFormat, << EI_Specification::val << ",\n" << EI_TokenizerMsg::val << "\n" );

	TYPEDEF_ERR_INFO(EI_ColumnIndex, unsigned int);
    DECLARE_INPUT_EXCEPTION(ExcMultipleDefinitionCsvColumn, << "Multiple definition of column with index " << EI_ColumnIndex::qval
    		<< " in included CSV file:\n" << EI_File::val << ",\n" );


    /**
     * @brief Read a storage from input stream.
     *
     * Parameter @p root_type provides input type tree declaration. See @p read_from_stream for details.
     */
    ReaderToStorage(const FilePath &in_file, Type::TypeBase &root_type);

    /// Read a storage from string (e.g. complex default value).
    ReaderToStorage( const string &default_str, Type::TypeBase &root_type, FileFormat format);

    /**
     * @brief Returns the root accessor.
     *
     * The template type \p T should correspond to the kind of the input type at root of the declaration tree.
     */
    template <class T>
    T get_root_interface() const;


    /**
     * @brief Default constructor.
     *
     * Provides common initialization for public constructors.
     */
    ReaderToStorage();

    /**
     * @brief This method actually reads the given stream \p in
     *
     * Checks the data just read against the declaration tree given by \p root_type, and
     * store the data into private storage tree using \p StorageBase classes.
     */
    void read_stream(istream &in, const Type::TypeBase &root_type, FileFormat format);

protected:

    /// Allow to sign reading the part of input file with atypical algorithm of reading.
    typedef enum {
    	none,         ///< Ordinary part of input.
		transposed,   ///< Processing transposed part of input tree.
		csv_include   ///< Processing part of input tree included in CSV file.
    } TryRead;

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

    /// Getter for root of the storage tree.
    StorageBase *get_storage();


    /// Storage of the read and checked input data
    StorageBase *storage_;

    /// Root of the declaration tree of the data in the storage.
    const Type::TypeBase *root_type_;

    friend class Type::Default;
    friend class ReaderInternalBase;

};











} // namespace Input



#endif /* READER_TO_STORAGE_HH_ */
