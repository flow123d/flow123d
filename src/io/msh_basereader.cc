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
 * @file    msh_basereader.cc
 * @brief
 * @author  dalibor
 */


#include "io/msh_basereader.hh"
#include "io/msh_gmshreader.h"
#include "io/msh_vtkreader.hh"
#include "system/sys_profiler.hh"


BaseMeshReader::BaseMeshReader(const Input::Record &mesh_rec)
: tok_(mesh_rec.val<FilePath>("mesh_file")) {
	current_cache_ = new ElementDataCache<double>();
	input_mesh_rec_ = mesh_rec;
}

BaseMeshReader::BaseMeshReader(std::istream &in)
: tok_(in) {
	current_cache_ = new ElementDataCache<double>();
}

std::shared_ptr< BaseMeshReader > BaseMeshReader::reader_factory(const Input::Record &mesh_rec) {
	FilePath file_path = mesh_rec.val<FilePath>("mesh_file");
	std::shared_ptr<BaseMeshReader> reader_ptr;
	try {
		if ( file_path.extension() == ".msh" ) {
			reader_ptr = std::make_shared<GmshMeshReader>(mesh_rec);
		} else if ( file_path.extension() == ".vtu" ) {
			reader_ptr = std::make_shared<VtkMeshReader>(mesh_rec);
		} else {
			THROW(ExcWrongExtension() << EI_FileExtension(file_path.extension()) << EI_MeshFile((string)file_path) );
		}
    } INPUT_CATCH(FilePath::ExcFileOpen, FilePath::EI_Address_String, mesh_rec)
	catch (ExceptionBase const &e) {
		throw;
	}
	return reader_ptr;

}

Mesh * BaseMeshReader::read_mesh() {
    START_TIMER("GMSHReader - read mesh");

	Input::Array region_list;
    Mesh * mesh = new Mesh( input_mesh_rec_ );
    this->read_physical_names(mesh);
	if (input_mesh_rec_.opt_val("regions", region_list)) {
		mesh->read_regions_from_input(region_list);
	}
    this->read_raw_mesh(mesh);
    mesh->setup_topology();
    mesh->check_and_finish();
    return mesh;
}

void BaseMeshReader::read_raw_mesh(Mesh* mesh) {
	ASSERT_PTR(mesh).error("Argument mesh is NULL.\n");
    tok_.set_position( Tokenizer::Position() );
    read_nodes(mesh);
    read_elements(mesh);
}

template<typename T>
typename ElementDataCache<T>::ComponentDataPtr BaseMeshReader::get_element_data( std::string field_name, double time,
		unsigned int n_entities, unsigned int n_components, bool &actual, std::vector<int> const & el_ids, unsigned int component_idx) {
	ASSERT(has_compatible_mesh_)
			.error("Vector of mapping VTK to GMSH element is not initialized. Did you call check_compatible_mesh?");

    MeshDataHeader actual_header = this->find_header(time, field_name);
    if ( !current_cache_->is_actual(actual_header.time, field_name) ) {
    	unsigned int size_of_cache; // count of vectors stored in cache

	    // check that the header is valid, try to correct
	    if (actual_header.n_entities != n_entities) {
	    	WarningOut().fmt("In file '{}', '{}' section for field '{}', time: {}.\nWrong number of entities: {}, using {} instead.\n",
	                tok_.f_name(), data_section_name_, field_name, actual_header.time, actual_header.n_entities, n_entities);
	        // actual_header.n_entities=n_entities;
	    }

	    if (n_components == 1) {
	    	// read for MultiField to 'n_comp' vectors
	    	// or for Field if ElementData contains only one value
	    	size_of_cache = actual_header.n_components;
	    }
	    else {
	    	// read for Field if more values is stored to one vector
	    	size_of_cache = 1;
	    	if (actual_header.n_components != n_components) {
	    		WarningOut().fmt("In file '{}', '{}' section for field '{}', time: {}.\nWrong number of components: {}, using {} instead.\n",
		                tok_.f_name(), data_section_name_, field_name, actual_header.time, actual_header.n_components, n_components);
		        actual_header.n_components=n_components;
	    	}
	    }

	    // set new cache
	    delete current_cache_;
	    current_cache_ = new ElementDataCache<T>(actual_header, size_of_cache, n_components*n_entities);

	    this->read_element_data(*current_cache_, actual_header, size_of_cache, n_components, el_ids);
	    actual = true; // use input header to indicate modification of @p data buffer
	}

    if (component_idx == std::numeric_limits<unsigned int>::max()) component_idx = 0;
	return static_cast< ElementDataCache<T> *>(current_cache_)->get_component_data(component_idx);
}


// explicit instantiation of template methods
#define MESH_READER_GET_ELEMENT_DATA(TYPE) \
template typename ElementDataCache<TYPE>::ComponentDataPtr BaseMeshReader::get_element_data<TYPE>(std::string field_name, double time, \
	unsigned int n_entities, unsigned int n_components, bool &actual, std::vector<int> const & el_ids, unsigned int component_idx);

MESH_READER_GET_ELEMENT_DATA(int);
MESH_READER_GET_ELEMENT_DATA(unsigned int);
MESH_READER_GET_ELEMENT_DATA(double);

