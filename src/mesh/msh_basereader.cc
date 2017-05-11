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


#include "mesh/msh_basereader.hh"


BaseMeshReader::BaseMeshReader(const FilePath &file_name)
: tok_(file_name) {
	current_cache_ = new ElementDataCache<double>();
}

BaseMeshReader::BaseMeshReader(std::istream &in)
: tok_(in) {
	current_cache_ = new ElementDataCache<double>();
}

template<typename T>
typename ElementDataCache<T>::ComponentDataPtr BaseMeshReader::get_element_data( std::string field_name, double time,
		unsigned int n_entities, unsigned int n_components, bool &actual, std::vector<int> const & el_ids, unsigned int component_idx) {
    MeshDataHeader actual_header = this->find_header(time, field_name);
    if ( !current_cache_->is_actual(actual_header.time, field_name) ) {
    	unsigned int size_of_cache; // count of vectors stored in cache

	    // check that the header is valid, try to correct
	    if (actual_header.n_entities != n_entities) {
	    	WarningOut().fmt("In file '{}', '$ElementData' section for field '{}', time: {}.\nWrong number of entities: {}, using {} instead.\n",
	                tok_.f_name(), field_name, actual_header.time, actual_header.n_entities, n_entities);
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
	    		WarningOut().fmt("In file '{}', '$ElementData' section for field '{}', time: {}.\nWrong number of components: {}, using {} instead.\n",
		                tok_.f_name(), field_name, actual_header.time, actual_header.n_components, n_components);
		        actual_header.n_components=n_components;
	    	}
	    }

	    // set new cache
	    delete current_cache_;
	    typename ElementDataCache<T>::CacheData data_cache = ElementDataCache<T>::create_data_cache(size_of_cache, n_components*n_entities);
	    current_cache_ = new ElementDataCache<T>(actual_header.time, actual_header.field_name, data_cache);

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

