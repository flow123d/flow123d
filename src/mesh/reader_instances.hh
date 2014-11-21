/*
 * reader_instances.hh
 *
 *  Created on: Feb 13, 2014
 *      Author: jb
 */

#ifndef READER_INSTANCES_HH_
#define READER_INSTANCES_HH_


#include "mesh/msh_gmshreader.h"
#include "system/file_path.hh"


/**
 * Auxiliary class to map filepaths to instances of readers.
 */
class ReaderInstances {
public:
	typedef std::map< string, std::shared_ptr<GmshMeshReader> > ReaderTable;

	/// Returns singleton instance
	static ReaderInstances * instance();

	/**
	 * Creates new reader if doesn't exist and returns its.
	 */
	std::shared_ptr<GmshMeshReader> create_reader(const FilePath &file_path);

	/**
	 * Returns mesh reader of get filepath. If reader doesn't exist, creates its.
	 */
	std::shared_ptr<GmshMeshReader> get_reader(const FilePath &file_path);

private:
	/// Constructor
	ReaderInstances() {};

	ReaderTable reader_table_;
};


#endif /* READER_INSTANCES_HH_ */
