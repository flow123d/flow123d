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
 * @file    msh_vtkreader.hh
 * @brief
 * @author  dalibor
 */

#ifndef MSH_VTK_READER_HH
#define	MSH_VTK_READER_HH

#include <string>
#include <istream>
#include <pugixml.hpp>

#include "system/file_path.hh"

class VtkMeshReader {
public:
    /**
     * Construct the VTK format reader from given filename.
     * This opens the file for reading.
     */
	VtkMeshReader(const FilePath &file_name);

    /**
     * Construct the VTK format reader from given input stream.
     * The input stream should be correctly opened. To get correct information about
     * line numbers there should be no previous reading from the stream.
     */
	VtkMeshReader(std::istream &in);

	/// Return count of nodes
	inline unsigned int n_nodes() const {
		return n_nodes_;
	}

	/// Return count of elements
	inline unsigned int n_elements() const {
		return n_elements_;
	}

protected:
	/// Empty constructor only for tests.
	VtkMeshReader() {}

	/// Set count of nodes and elements.
	void read_nodes_elms_count();

	pugi::xml_document doc_;
    pugi::xml_parse_result parse_result_;

    /// count of nodes
    unsigned int n_nodes_;

    /// count of elements
    unsigned int n_elements_;

};

#endif	/* MSH_VTK_READER_HH */

