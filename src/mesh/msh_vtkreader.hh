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
	/// Possible data sections in UnstructuredGrid - Piece node.
	enum DataSections {
	    points, cells, cell_data
	};

	/// Type of data formats - ascii or appended.
	enum DataFormat {
	    ascii, appended
	};

	/// Type of VTK data (only supported formats!)
	enum DataType {
	    uint32, uint64, float64
	};

	/// Attributes of DataArray section
	struct DataArrayAttributes {
		std::string name_;          ///< Name of DataArray
		DataType type_;             ///< Type (only UInt32, UInt64 and Float64 are supported)
		unsigned int n_components_; ///< NumberOfComponents (default value is 1)
		DataFormat format_;         ///< Format (ascii or appended)
		unsigned int offset_;       ///< Offset of data (only for appended format)
	};

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

	/**
	 * Get XML node in UnstructuredGrid part of VTK file.
	 * @param data_section    Section where node is located.
	 * @param data_array_name Attribute "Name" of DataArray tag (not used for Point section)
	 */
	DataArrayAttributes get_data_array_attr(DataSections data_section, std::string data_array_name = "");

	/// Return count of nodes
	inline unsigned int n_nodes() const {
		return n_nodes_;
	}

	/// Return count of elements
	inline unsigned int n_elements() const {
		return n_elements_;
	}

protected:
	/// Get DataType by value of string
	DataType get_data_type(std::string type_str, bool only_integral = false);

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

    /// header type of VTK file
    std::string header_type_;

    /// compressed data of AppendedData section
    bool compressed_;

};

#endif	/* MSH_VTK_READER_HH */

