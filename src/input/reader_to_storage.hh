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
    DECLARE_INPUT_EXCEPTION( ExcForbiddenAbstractTag, << "Forbidden tag " << EI_Tag::qval << " of Abstract type.\n" );


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
    	type_int, type_double, type_bool, type_string
    } IncludeDataTypes;

    /// Data of one column of including CSV file.
    struct IncludeCsvData {
    	IncludeDataTypes data_type;
    	vector<unsigned int> storage_indexes;
    };

    /// Getter for root of the storage tree.
    StorageBase *get_storage();


    /**
     * @brief Create storage of given @p type.
     *
     * Check correctness of the input given by json_spirit or YAML-cpp node at head() of PathBase @p p
     * against type specification @p type. Die on input error (and return NULL).
     * For correct input, creates the storage tree and returns pointer to its root node.
     */
    StorageBase * make_storage(PathBase &p, const Type::TypeBase *type);

    StorageBase * make_storage(PathBase &p, const Type::Record *record);         ///< Create storage of Type::Record type
    StorageBase * make_storage(PathBase &p, const Type::Abstract *abstr_rec);    ///< Create storage of Type::Abstract type
    StorageBase * make_storage(PathBase &p, const Type::Array *array);           ///< Create storage of Type::Array type
    StorageBase * make_storage(PathBase &p, const Type::Tuple *tuple);		 ///< Create storage of Type::Tuple type
    StorageBase * make_storage(PathBase &p, const Type::Selection *selection);   ///< Create storage of Type::Selection type
    StorageBase * make_storage(PathBase &p, const Type::Bool *bool_type);        ///< Create storage of Type::Bool type
    StorageBase * make_storage(PathBase &p, const Type::Integer *int_type);      ///< Create storage of Type::Integer type
    StorageBase * make_storage(PathBase &p, const Type::Double *double_type);    ///< Create storage of Type::Double type
    StorageBase * make_storage(PathBase &p, const Type::String *string_type);    ///< Create storage of Type::String type

    /// Apply transposition and create storage of Type::Array type
    StorageBase * make_transposed_storage(PathBase &p, const Type::TypeBase *type);

    /// Apply conversion to one element storage of Type::Array type
    StorageBase * make_autoconversion_array_storage(PathBase &p, const Type::Array *array, StorageBase *item);

    /// Apply automatic conversion of Type::Record type
    StorageBase * record_automatic_conversion(PathBase &p, const Type::Record *record);

    /// Apply automatic conversion of Type::Abstract type
    StorageBase * abstract_automatic_conversion(PathBase &p, const Type::Abstract *abstr_rec);

    /// Dispatch according to @p type and create corresponding storage from the given string.
    StorageBase * make_storage_from_default( const string &dflt_str, std::shared_ptr<Type::TypeBase> type);

    /// Create storage of included YAML or JSON input file
    StorageBase * make_include_storage(PathBase &p, const Type::Record *record);

    /// Create storage of included CSV input file
    StorageBase * make_include_csv_storage(PathBase &p, const Type::Array *array);

    /// Helper method. Get string value of included file or throw exception if reading failed.
    std::string get_included_file(PathBase &p);

    /// Set storage of simple input type with value given from CSV file.
    void set_storage_from_csv(unsigned int column_index, StorageBase * item_storage, StorageBase * new_storage);


    /// Storage of the read and checked input data
    StorageBase *storage_;

    /// Root of the declaration tree of the data in the storage.
    const Type::TypeBase *root_type_;

    /**
     * @brief Flag signed that "special" part of input tree is processed.
     *
     * A. <b>Transposed part of input tree:</b>
     *
     * We set this flag if input tree contains another type at position where Array
     * is expected. This type must correspond with type_of_value of Array.
     *
     * Subsequently:
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
     *
     *
     * B. <b>Include of CSV file to input tree:</b>
     *
     * Not implemented yet.
     */
    TryRead try_read_;

    /// Index of processed item in transposed part of input tree.
    unsigned int transpose_index_;

    /// Helper vector what allows check sizes of all transposed Arrays.
    vector<unsigned int> transpose_array_sizes_;

    /// Helper vector which contains actual indexes of subtree imported in CSV file.
    vector<unsigned int> csv_storage_indexes_;

    /// Map of columns in CSV file to storage of subtree
    map<unsigned int, IncludeCsvData> csv_columns_map_;

    friend class Type::Default;

};











} // namespace Input



#endif /* READER_TO_STORAGE_HH_ */
