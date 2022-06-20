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
#include "system/index_types.hh" // for LongIdx

class BaseMeshReader;
class Mesh;
class EquivalentMeshMap;



/**
 * Auxiliary class to map filepaths to instances of readers.
 */
class ReaderCache {
public:
	struct ReaderData {
		/// Constructor
		ReaderData() : target_mesh_element_map_(nullptr) {};

		std::shared_ptr<BaseMeshReader> reader_;
		std::shared_ptr<Mesh> mesh_;
		std::shared_ptr<EquivalentMeshMap> target_mesh_element_map_;
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
	 * Fill element id vectors of reader without checking compatibility.
	 */
	static void get_element_ids(const FilePath &file_path, const Mesh &mesh);

	/**
	 * Returns shared vector mapping elements to target mesh.
	 *
	 * Reader and appropriate input data mesh are given by FilePath.
	 * If map is not created method check_compatible_mesh of \p computational_mesh is called.
	 */
    static std::shared_ptr<EquivalentMeshMap> eqivalent_mesh_map(const FilePath &file_path,
                                                                          Mesh *computational_mesh);

    static std::shared_ptr<EquivalentMeshMap> identic_mesh_map(const FilePath &file_path,
                                                                          Mesh *computational_mesh);

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
