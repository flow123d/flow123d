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
 * @file    msh_pvdreader.cc
 * @brief
 * @author  dalibor
 */


#include <iostream>
#include <vector>
#include <pugixml.hpp>

#include "msh_pvdreader.hh"
#include "msh_vtkreader.hh"



PvdMeshReader::PvdMeshReader(const FilePath &file_name)
: BaseMeshReader(file_name)
{
    data_section_name_ = "DataArray";
    has_compatible_mesh_ = false;
    pvd_path_dir_ = file_name.parent_path();
	make_header_table();
}


PvdMeshReader::~PvdMeshReader()
{
	for (auto file_data : file_list_) {
		if (file_data.reader != nullptr) delete file_data.reader;
	}
}


void PvdMeshReader::read_physical_names(Mesh * mesh) {
	// will be implemented later
}


void PvdMeshReader::read_nodes(Mesh * mesh) {
	ASSERT(false).error("Reading of VTK mesh is not supported yet!");
	// will be implemented later
}


void PvdMeshReader::read_elements(Mesh * mesh) {
	// will be implemented later
}


void PvdMeshReader::read_element_data(ElementDataCacheBase &data_cache, MeshDataHeader actual_header, unsigned int n_components,
		bool boundary_domain) {

	ASSERT(!boundary_domain).error("Reading PVD data of boundary elements is not supported yet!\n");

}


void PvdMeshReader::check_compatible_mesh(Mesh &mesh) {
	// will be implemented later
}


void PvdMeshReader::make_header_table() {
	pugi::xml_document doc;;
	doc.load_file( tok_.f_name().c_str() );
	pugi::xml_node node = doc.child("VTKFile").child("Collection");

	double last_time = -numeric_limits<double>::infinity(); // check ascending order of times
	std::vector<std::string> sub_paths; //allow construct paths of VTK files
	sub_paths.resize(2);
	sub_paths[0] = pvd_path_dir_;

	// read PVD data
	for (pugi::xml_node subnode = node.child("DataSet"); subnode; subnode = subnode.next_sibling("DataSet")) {
		double time = subnode.attribute("timestep").as_double();
		if (time <= last_time) {
			WarningOut().fmt("Wrong time order in PVD file '{}', time '{}'. Skipping this time step.\n", tok_.f_name(), time );
		} else {
			sub_paths[1] = subnode.attribute("file").as_string();
			FilePath vtk_path(sub_paths, FilePath::input_file);
			last_time = time;
			file_list_.push_back( VtkFileData(time, vtk_path) );
		}
	}
}


MeshDataHeader & PvdMeshReader::find_header(double time, std::string field_name) {
	auto comp = [](double t, const VtkFileData &a) {
		return t * (1.0 + 2.0*numeric_limits<double>::epsilon()) < a.time;
	};

	// find iterator to data of VTK file
	std::vector<VtkFileData>::iterator list_it = std::upper_bound(file_list_.begin(),
			file_list_.end(),
			time,
			comp);
	--list_it;

	// check if VTK reader exists and eventually creates its
	if (!list_it->reader) {
		list_it->reader = new VtkMeshReader(list_it->file_name, this->element_data_values_, list_it->time);
	}

	return list_it->reader->find_header(time, field_name);
}

