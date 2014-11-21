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

std::shared_ptr<GmshMeshReader> ReaderInstances::create_reader(const FilePath &file_path) {
	ReaderTable::iterator it = reader_table_.find( string(file_path) );
	if (it == reader_table_.end()) {
		std::shared_ptr<GmshMeshReader> reader_ptr = std::make_shared<GmshMeshReader>(file_path);
		reader_table_.insert( std::pair<string, std::shared_ptr<GmshMeshReader>>(string(file_path), reader_ptr) );
		return reader_ptr;
	} else {
		xprintf(Warn, "Instance of reader %s already exists, it can't be created!\n", string(file_path).c_str());
		return (*it).second;
	}
}

std::shared_ptr<GmshMeshReader> ReaderInstances::get_reader(const FilePath &file_path) {
	ReaderTable::iterator it = reader_table_.find( string(file_path) );
	if (it == reader_table_.end()) {
		xprintf(Warn, "Instance of reader %s doesn't exist, it will be created.\n", string(file_path).c_str());
		std::shared_ptr<GmshMeshReader> reader_ptr = std::make_shared<GmshMeshReader>(file_path);
		reader_table_.insert( std::pair<string, std::shared_ptr<GmshMeshReader>>(string(file_path), reader_ptr) );
		return reader_ptr;
	} else {
		return (*it).second;
	}
}
