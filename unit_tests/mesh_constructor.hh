/*
 * mesh_constructor.hh
 *
 *  Created on: Jul 4, 2016
 *      Author: jb
 */

#ifndef UNIT_TESTS_MESH_CONSTRUCTOR_HH_
#define UNIT_TESTS_MESH_CONSTRUCTOR_HH_

#include <iostream>

#include "input/reader_to_storage.hh"
#include "input/type_record.hh"
#include "mesh/mesh.h"


namespace IT = Input::Type;

/**
 * Construct mesh.
 */
Mesh * mesh_constructor(const std::string &input_str="{mesh_file=\"\"}",
		Input::FileFormat format = Input::FileFormat::format_JSON, MPI_Comm com = MPI_COMM_WORLD) {

	istringstream is(input_str);
    Input::ReaderToStorage reader;
    IT::Record &in_rec = const_cast<IT::Record &>(Mesh::get_input_type());
    in_rec.finish();
    reader.read_stream(is, in_rec, format);

    return new Mesh(reader.get_root_interface<Input::Record>(), com);
}


#endif /* UNIT_TESTS_MESH_CONSTRUCTOR_HH_ */
