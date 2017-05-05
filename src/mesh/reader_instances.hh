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
 * @file    reader_instances.hh
 * @brief   
 */

#ifndef READER_INSTANCES_HH_
#define READER_INSTANCES_HH_


#include "mesh/msh_basereader.hh"
#include "system/file_path.hh"


/**
 * Auxiliary class to map filepaths to instances of readers.
 */
class ReaderInstances {
public:
	TYPEDEF_ERR_INFO(EI_MeshFile, std::string);
	TYPEDEF_ERR_INFO(EI_FileExtension, std::string);
	DECLARE_EXCEPTION(ExcWrongExtension,
			<< "Unsupported extension " << EI_FileExtension::qval << " of the input file: " << EI_MeshFile::qval);

	typedef std::map< string, std::shared_ptr<BaseMeshReader> > ReaderTable;

	/// Returns singleton instance
	static ReaderInstances * instance();

	/**
	 * Returns mesh reader of get filepath. If reader doesn't exist, creates its.
	 */
	std::shared_ptr<BaseMeshReader> get_reader(const FilePath &file_path);

private:
	/// Constructor
	ReaderInstances() {};

	/// Table of readers
	ReaderTable reader_table_;
};


#endif /* READER_INSTANCES_HH_ */
