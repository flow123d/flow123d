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

#include "input/input_type.hh"

#include "input/input_exception.hh"
#include "input/storage.hh"
#include "input/path_base.hh"


namespace Input {



/// Possible formats of input files.
typedef enum  {
    format_JSON,
    format_YAML
} FileFormat;



/// Enum of possible input types. Values in @p json_type_names must be stored in same order.
typedef enum {
	obj_type, array_type, str_type, bool_type, int_type, real_type, null_type, scalar_type, undef_type
} ValueTypes;



/**
 *  @brief Reader for (slightly) modified JSON files.
 *
 *  This class implements a reader of modified JSON file format. The modifications include
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


    /**
     * Read a storage from input stream. Parameter @p root_type
     * provides input type tree declaration. See @p read_from_stream for details.
     */
    ReaderToStorage(const FilePath &in_file, const Type::TypeBase &root_type);

    /**
     * Read a storage from string (e.g. complex default value).
     */
    ReaderToStorage( const string &default_str, const Type::TypeBase &root_type, FileFormat format);

    /**
     * Returns the root accessor. The template type \p T should correspond
     * to the kind of the input type at root of the declaration tree.
     */
    template <class T>
    T get_root_interface() const;


protected:

    /**
     * Default constructor.
     * Provides common initialization for public constructors.
     */
    ReaderToStorage();

    /**
     * This method actually reads the given stream \p in, checks the data just read against the declaration tree given by \p root_type, and
     * store the data into private storage tree using \p StorageBase classes.
     */
    void read_stream(istream &in, const Type::TypeBase &root_type, FileFormat format);

    /**
     * Getter for root of the storage tree.
     */
    StorageBase *get_storage()
    { return storage_;}


    /**
     * Check correctness of the input given by json_spirit node at head() of PathJSON @p p
     * against type specification @p type. Die on input error (and return NULL).
     * For correct input, creates the storage tree and returns pointer to its root node.
     */
    StorageBase * make_storage(PathBase &p, const Type::TypeBase *type);

    StorageBase * make_storage(PathBase &p, const Type::Record *record);
    StorageBase * make_storage(PathBase &p, const Type::AbstractRecord *abstr_rec);
    StorageBase * make_storage(PathBase &p, const Type::Array *array);

    StorageBase * make_selection_storage_without_catch(PathBase &p, const Type::Selection *selection);
    StorageBase * make_storage(PathBase &p, const Type::Selection *selection);
    StorageBase * make_storage(PathBase &p, const Type::Bool *bool_type);
    StorageBase * make_storage(PathBase &p, const Type::Integer *int_type);
    StorageBase * make_storage(PathBase &p, const Type::Double *double_type);
    StorageBase * make_storage(PathBase &p, const Type::String *string_type);

    StorageBase * record_automatic_conversion(PathBase &p, const Type::Record *record);
    StorageBase * abstract_rec_automatic_conversion(PathBase &p, const Type::AbstractRecord *abstr_rec);

    /**
     * Dispatch according to @p type and create corresponding storage from the given string.
     */
    StorageBase * make_storage_from_default( const string &dflt_str, boost::shared_ptr<Type::TypeBase> type);


    /// Storage of the read and checked input data
    StorageBase *storage_;

    /// Root of the declaration tree of the data in the storage.
    const Type::TypeBase *root_type_;

    friend class Type::Default;

};











} // namespace Input



#endif /* READER_TO_STORAGE_HH_ */
