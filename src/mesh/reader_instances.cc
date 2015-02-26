/*
 * reader_instances.hh
 *
 *  Created on: Feb 13, 2014
 *      Author: jb
 */

#include "mesh/reader_instances.hh"

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
