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


/**
 * Initialize mesh by data structures represent mesh stored in UNIT_TESTS_SRC_DIR+/mesh/simplest_cube.msh
 */
void init_simplest_cube_mesh(Mesh *mesh) {
	PhysicalNamesDataTable physical_names_data = {
			{1, 37, "1D diagonal"},
			{2, 38, "2D XY diagonal"},
			{2, 101, ".top side"},
			{2, 102, ".bottom side"},
			{3, 39, "3D back"},
			{3, 40, "3D front"}
	};
	NodeDataTable node_data = {
			{1, {1, 1, 1} },
			{2, {-1, 1, 1} },
			{3, {-1, -1, 1} },
			{4, {1, -1, 1} },
			{5, {1, -1, -1} },
			{6, {-1, -1, -1} },
			{7, {1, 1, -1} },
			{8, {-1, 1, -1} }
	};
	ElementDataTable element_data = {
			{ 1, 1, 37, 0, 7, 3 },
			{ 2, 2, 38, 0, 6, 3, 7 },
			{ 3, 2, 38, 0, 3, 1, 7 },
			{ 4, 3, 39, 0, 3, 7, 1, 2 },
			{ 5, 3, 39, 0, 3, 7, 2, 8 },
			{ 6, 3, 39, 0, 3, 7, 8, 6 },
			{ 7, 3, 40, 0, 3, 7, 6, 5 },
			{ 8, 3, 40, 0, 3, 7, 5, 4 },
			{ 9, 3, 40, 0, 3, 7, 4, 1 },
			{10, 2, 101, 0, 1, 2, 3 },
			{11, 2, 101, 0, 1, 3, 4 },
			{12, 2, 102, 0, 6, 7, 8 },
			{13, 2, 102, 0, 7, 6, 5 }
	};
	mesh->init_from_input(physical_names_data, node_data, element_data);
}


#endif /* UNIT_TESTS_MESH_CONSTRUCTOR_HH_ */
