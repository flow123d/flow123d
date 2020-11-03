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
 * @file    reader_cache.hh
 * @brief   
 */

#ifndef READER_CACHE_HH_
#define READER_CACHE_HH_


#include <map>                  // for map, map<>::value_compare
#include <memory>               // for shared_ptr
#include <string>               // for string
#include "system/file_path.hh"  // for FilePath

class BaseMeshReader;
class Mesh;



/**
 * Auxiliary class to map filepaths to instances of readers.
 */
class ReaderCache {
public:
	struct ReaderData {
		std::shared_ptr<BaseMeshReader> reader_;
		std::shared_ptr<Mesh> mesh_;
	};

	typedef std::map< string, ReaderData > ReaderTable;

	/**
	 * Returns reader of given FilePath.
	 */
	static std::shared_ptr<BaseMeshReader> get_reader(const FilePath &file_path);

	/**
	 * Returns mesh of given FilePath.
	 */
	static std::shared_ptr<Mesh> get_mesh(const FilePath &file_path);

	/**
	 * Check if nodes and elements of reader mesh are compatible with \p mesh and fill element id vectors of reader.
	 *
	 * OBSOLETE method - will be removed.
	 */
	static bool check_compatible_mesh(const FilePath &file_path, Mesh &mesh);

	/**
	 * Fill element id vectors of reader without checking compatibility.
	 *
	 * OBSOLETE method - will be removed or change.
	 */
	static void get_element_ids(const FilePath &file_path, const Mesh &mesh);

private:
	/// Returns singleton instance
	static ReaderCache * instance();

	/// Constructor
	ReaderCache() {};

	/// Returns instance of given FilePath. If reader doesn't exist, creates new ReaderData object.
	static ReaderTable::iterator get_reader_data(const FilePath &file_path);

	/// Table of readers
	ReaderTable reader_table_;
};


#endif /* READER_CACHE_HH_ */
