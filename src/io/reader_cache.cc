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
 * @file    reader_cache.cc
 * @brief   
 */

#include "io/reader_cache.hh"
#include "io/msh_basereader.hh"
#include "io/msh_gmshreader.h"
#include "io/msh_vtkreader.hh"
#include "io/msh_pvdreader.hh"
#include "mesh/mesh.h"
#include "input/accessors.hh"


/***********************************************************************************************
 * Implementation of ReaderCache
 */

ReaderCache * ReaderCache::instance() {
	static ReaderCache *instance = new ReaderCache;
	return instance;
}

std::shared_ptr<BaseMeshReader> ReaderCache::get_reader(const FilePath &file_path) {
	return ReaderCache::get_reader_data(file_path)->second.reader_;
}

std::shared_ptr<Mesh> ReaderCache::get_mesh(const FilePath &file_path) {
	auto it = ReaderCache::get_reader_data(file_path);
	// Create and fill mesh if doesn't exist
	if ( (*it).second.mesh_ == nullptr ) {
		(*it).second.mesh_ = std::make_shared<Mesh>( Input::Record() );
		(*it).second.reader_->read_physical_names( (*it).second.mesh_.get() );
		(*it).second.reader_->read_raw_mesh( (*it).second.mesh_.get() );
		//(*it).second.reader_->check_compatible_mesh( *((*it).second.mesh_) );

	}
	return (*it).second.mesh_;
}

ReaderCache::ReaderTable::iterator ReaderCache::get_reader_data(const FilePath &file_path) {
	ReaderTable::iterator it = ReaderCache::instance()->reader_table_.find( string(file_path) );
	if (it == ReaderCache::instance()->reader_table_.end()) {
		ReaderData reader_data;
		if ( file_path.extension() == ".msh" ) {
			reader_data.reader_ = std::make_shared<GmshMeshReader>(file_path);
		} else if ( file_path.extension() == ".vtu" ) {
			reader_data.reader_ = std::make_shared<VtkMeshReader>(file_path);
		} else if ( file_path.extension() == ".pvd" ) {
			reader_data.reader_ = std::make_shared<PvdMeshReader>(file_path);
		} else {
			THROW(BaseMeshReader::ExcWrongExtension()
				<< BaseMeshReader::EI_FileExtension(file_path.extension()) << BaseMeshReader::EI_MeshFile((string)file_path) );
		}
		ReaderCache::instance()->reader_table_.insert( std::pair<string, ReaderData>(string(file_path), reader_data) );
		it = ReaderCache::instance()->reader_table_.find( string(file_path) );
	}
	return it;
}

bool ReaderCache::check_compatible_mesh(const FilePath &file_path, Mesh &mesh) {
	auto mesh_ptr = ReaderCache::get_mesh(file_path);
	auto reader_ptr = ReaderCache::get_reader(file_path);
	reader_ptr->has_compatible_mesh_ = true;
	return mesh_ptr->check_compatible_mesh(mesh, reader_ptr->bulk_elements_id_, reader_ptr->boundary_elements_id_);
}

void ReaderCache::get_element_ids(const FilePath &file_path, const Mesh &mesh) {
	auto reader_ptr = ReaderCache::get_reader(file_path);
	reader_ptr->has_compatible_mesh_ = true;
	mesh.elements_id_maps(reader_ptr->bulk_elements_id_, reader_ptr->boundary_elements_id_);
}
