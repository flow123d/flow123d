/*
 * reader_to_storage_constructor.hh
 *
 *  Created on: Jul 4, 2016
 *      Author: jb
 */

#ifndef UNIT_TESTS_READER_ROOT_INTERFACE_HH_
#define UNIT_TESTS_READER_ROOT_INTERFACE_HH_

#include <iostream>

#include "input/reader_to_storage.hh"
#include "input/type_base.hh"


namespace IT = Input::Type;

/**
 * Get root interface of input tree represented by root type.
 */
template <class T>
T reader_root_interface(const std::string &input_str, IT::TypeBase &root_type,
		Input::FileFormat format = Input::FileFormat::format_JSON) {

	istringstream is(input_str);
	root_type.finish();
    Input::ReaderToStorage reader;
    reader.read_stream(is, root_type, format);

    return reader.get_root_interface<T>();;
}


#endif /* UNIT_TESTS_READER_ROOT_INTERFACE_HH_ */
