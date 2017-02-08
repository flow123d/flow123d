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

#include "mesh/reader_instances.hh"
#include "input/accessors.hh"


ReaderInstance::ReaderData ReaderInstance::get_instance(const FilePath &file_path) {
	static ReaderInstance *instance = new ReaderInstance;

	ReaderTable::iterator it = instance->reader_table_.find( string(file_path) );
	if (it == instance->reader_table_.end()) {
		ReaderData reader_data;
		reader_data.reader_ = std::make_shared<GmshMeshReader>(file_path);
		reader_data.mesh_ = std::make_shared<Mesh>( Input::Record() );
		reader_data.reader_->read_mesh( reader_data.mesh_.get() );
		instance->reader_table_.insert( std::pair<string, ReaderData>(string(file_path), reader_data) );
		return reader_data;
	} else {
		return (*it).second;
	}
}


ReaderInstances * ReaderInstances::instance() {
	static ReaderInstances *instance = new ReaderInstances;
	return instance;
}

std::shared_ptr<GmshMeshReader> ReaderInstances::get_reader(const FilePath &file_path) {
	ReaderTable::iterator it = reader_table_.find( string(file_path) );
	if (it == reader_table_.end()) {
		std::shared_ptr<GmshMeshReader> reader_ptr = std::make_shared<GmshMeshReader>(file_path);
		reader_table_.insert( std::pair<string, std::shared_ptr<GmshMeshReader>>(string(file_path), reader_ptr) );
		return reader_ptr;
	} else {
		return (*it).second;
	}
}
