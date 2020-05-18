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
 * @file    reader_to_storage.cc
 * @brief   
 */

#include <cstdint>
#include <limits>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#include "reader_to_storage.hh"
#include "input/path_json.hh"
#include "input/path_yaml.hh"
#include "input/input_type.hh"
#include "input/accessors.hh"
#include "input/reader_internal.hh"
#include <stddef.h>                                    // for NULL
#include <boost/exception/detail/error_info_impl.hpp>  // for error_info
#include <ostream>                                     // for operator<<
#include <set>                                         // for set, _Rb_tree_...
#include <typeinfo>                                    // for type_info
#include <utility>                                     // for pair
#include <vector>                                      // for vector
#include "system/asserts.hh"                           // for Assert, ASSERT
#include "system/file_path.hh"                         // for FilePath, File...
#include "system/logger.hh"                            // for operator<<


namespace Input {
using namespace std;
using namespace internal;



/********************************************
 * Implementation of public part of ReaderToStorage
 */

ReaderToStorage::ReaderToStorage()
: storage_(nullptr),
  root_type_(nullptr)
{}



ReaderToStorage::ReaderToStorage(const FilePath &in_file, Type::TypeBase &root_type)
: ReaderToStorage()
{
	std::string fname = in_file;
	std::string extension = fname.substr(fname.find_last_of(".") + 1);
	FileFormat format;
	if (extension == "con") {
		format = FileFormat::format_JSON;
	} else if (extension == "yaml") {
		format = FileFormat::format_YAML;
	} else {
		THROW(ExcInputMessage() << EI_Message("Invalid extension of file " + fname + ".\nMust be 'con' or 'yaml'."));
	}

	try {
		std::ifstream in;
		in_file.open_stream(in);

		// finish of root_type ensures finish of whole IST
		root_type.finish();

		read_stream(in, root_type, format);
    } catch (Input::Exception &e ) {
      e << Input::ReaderInternalBase::EI_File(in_file); throw;
    }
}



ReaderToStorage::ReaderToStorage( const string &str, Type::TypeBase &root_type, FileFormat format)
: ReaderToStorage()
{
    // finish of root_type ensures finish of whole IST
    root_type.finish();

	try {
		istringstream is(str);
		read_stream(is, root_type, format);
	} catch (ReaderInternalBase::ExcNotJSONFormat &e) {
		e << ReaderInternalBase::EI_File("STRING: "+str); throw;
	}
}



StorageBase *ReaderToStorage::get_storage()
{
	return storage_;
}



void ReaderToStorage::read_stream(istream &in, const Type::TypeBase &root_type, FileFormat format)
{
	ASSERT(storage_==nullptr).error();

    PathBase * root_path;
	if (format == FileFormat::format_JSON) {
		root_path = new PathJSON(in);
	} else {
		root_path = new PathYAML(in);
	}

	// guarantee to delete root_path on function return even on exception
	std::unique_ptr<PathBase> root_path_ptr(root_path);

    root_type_ = &root_type;
	try {
		ReaderInternal ri;
	    storage_ = ri.read_storage(*root_path_ptr, root_type_);
	} catch (ReaderInternalBase::ExcInputError &e) {
		if (format == FileFormat::format_JSON) {
			e << ReaderInternalBase::EI_Format("JSON");
		} else {
			e << ReaderInternalBase::EI_Format("YAML");
		}
		throw;
	}

	ASSERT_PTR(storage_).error();
}






/********************************************88
 * Implementation
 */

template <class T>
T ReaderToStorage::get_root_interface() const
{
	ASSERT_PTR(storage_).error();

    Address addr(storage_, root_type_);
    // try to create an iterator just to check type
    Iterator<T>( *root_type_, addr, 0);

    auto tmp_root_type = static_cast<const typename T::InputType &>(*root_type_);
    return T( addr, tmp_root_type );
}


template ::Input::Record ReaderToStorage::get_root_interface<::Input::Record>() const;
template ::Input::Array ReaderToStorage::get_root_interface<::Input::Array>() const;
template ::Input::AbstractRecord ReaderToStorage::get_root_interface<::Input::AbstractRecord>() const;
template ::Input::Tuple ReaderToStorage::get_root_interface<::Input::Tuple>() const;
//template ReaderToStorage::get_root_interface<::Input::>()->::Input::AbstractRecord const;


} // namespace Input
