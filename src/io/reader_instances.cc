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

ReaderInstances * ReaderInstances::instance() {
	static ReaderInstances *instance = new ReaderInstances;
	return instance;
}

std::shared_ptr<BaseMeshReader> ReaderInstances::get_reader(const FilePath &file_name) {
	ReaderTable::iterator it = reader_table_.find( string(file_name) );
	if (it == reader_table_.end()) {
		std::shared_ptr<BaseMeshReader> reader_ptr;
		if ( file_name.extension() == ".msh" ) {
			reader_ptr = std::make_shared<GmshMeshReader>(file_name);
		} else if ( file_name.extension() == ".vtu" ) {
			reader_ptr = std::make_shared<VtkMeshReader>(file_name);
		} else {
			THROW(BaseMeshReader::ExcWrongExtension()
				<< BaseMeshReader::EI_FileExtension(file_name.extension()) << BaseMeshReader::EI_MeshFile((string)file_name) );
		}
		reader_table_.insert( std::pair<string, std::shared_ptr<BaseMeshReader>>(string(file_name), reader_ptr) );
		return reader_ptr;
	} else {
		return (*it).second;
	}
}
