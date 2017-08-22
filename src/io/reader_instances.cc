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
 * @file    reader_instances.cc
 * @brief   
 */

#include "io/reader_instances.hh"
#include "io/msh_gmshreader.h"
#include "io/msh_vtkreader.hh"
#include "io/msh_pvdreader.hh"
#include "input/accessors.hh"


/***********************************************************************************************
 * Implementation of ReaderInstance
 */

ReaderInstance * ReaderInstance::instance() {
	static ReaderInstance *instance = new ReaderInstance;
	return instance;
}

std::shared_ptr<BaseMeshReader> ReaderInstance::get_reader(const FilePath &file_path) {
	return ReaderInstance::get_instance(file_path).reader_;
}

std::shared_ptr<Mesh> ReaderInstance::get_mesh(const FilePath &file_path) {
	// First, we must check if reader instance item exists in table
	ReaderTable::iterator it = ReaderInstance::instance()->reader_table_.find( string(file_path) );
	if (it == ReaderInstance::instance()->reader_table_.end()) {
		ReaderInstance::get_instance(file_path); // add new reader instance to table
		it = ReaderInstance::instance()->reader_table_.find( string(file_path) );
	}
	// Create and fill mesh if doesn't exist
	if ( (*it).second.mesh_ == nullptr ) {
		(*it).second.mesh_ = std::make_shared<Mesh>( Input::Record() );
		(*it).second.reader_->read_raw_mesh( (*it).second.mesh_.get() );
		(*it).second.reader_->check_compatible_mesh( *((*it).second.mesh_) );

	}
	return (*it).second.mesh_;
}

ReaderInstance::ReaderData ReaderInstance::get_instance(const FilePath &file_path) {
	ReaderTable::iterator it = ReaderInstance::instance()->reader_table_.find( string(file_path) );
	if (it == ReaderInstance::instance()->reader_table_.end()) {
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
		ReaderInstance::instance()->reader_table_.insert( std::pair<string, ReaderData>(string(file_path), reader_data) );
		return reader_data;
	} else {
		return (*it).second;
	}
}

