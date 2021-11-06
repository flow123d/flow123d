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


#include "system/file_path.hh"
#include <sys/types.h>                 // for int64_t

#include <memory>                      // for shared_ptr
#include <string>                      // for string
#include "system/exceptions.hh"        // for operator<<, ExcStream, EI, TYP...
namespace Input { class PathBase; }
namespace Input { class StorageBase; }
namespace Input { namespace Type { class Abstract; } }
namespace Input { namespace Type { class Array; } }
namespace Input { namespace Type { class Bool; } }
namespace Input { namespace Type { class Double; } }
namespace Input { namespace Type { class Integer; } }
namespace Input { namespace Type { class Record; } }
namespace Input { namespace Type { class Selection; } }
namespace Input { namespace Type { class String; } }
namespace Input { namespace Type { class Tuple; } }
namespace Input { namespace Type { class TypeBase; } }

namespace Input { namespace Type { class Default; } }
namespace Input { namespace Type { class TypeBase; } }



namespace Input {

using namespace std;

//class ReaderInternalBase;
//class ReaderInternalCsvInclude;



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
