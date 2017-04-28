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

#include "mesh/msh_basereader.hh"
#include "system/file_path.hh"

class VtkMeshReader : public BaseMeshReader {
public:
	/// Possible data sections in UnstructuredGrid - Piece node.
	enum DataSections {
	    points, cells, cell_data
	};

	/// Type of data formats - ascii, binary or compressed with zLib.
	enum DataFormat {
	    ascii, binary_uncompressed, binary_zlib
	};

	/// Type of VTK data (value 'undefined' for empty value)
	enum DataType {
	    int8, uint8, int16, uint16, int32, uint32, int64, uint64, float32, float64, undefined
	};

	/// Attributes of DataArray section
	struct DataArrayAttributes {
		std::string field_name;     ///< Name of DataArray
		DataType type_;             ///< Type of data
		unsigned int n_components;  ///< NumberOfComponents (default value is 1)
		unsigned int offset_;       ///< Offset of data (only for appended format)
		std::string tag_value_;     ///< String value of tag (only for ascii format)
	};

	/**
	 * Map of DataArray sections in VTK file.
	 *
	 * For each field_name contains DataArrayAttributes.
	 */
	typedef typename std::map< std::string, DataArrayAttributes > HeaderTable;

	/**
     * Construct the VTK format reader from given filename.
     * This opens the file for reading.
     */
	VtkMeshReader(const FilePath &file_name);

	/// Destructor
	~VtkMeshReader();

    /**
	 * Find header of DataArray section of VTK file given by field_name.
	 *
	 * Note: \p time has no effect (it is only for continuity with GMSH reader).
	 */
	DataArrayAttributes find_header(double time, std::string field_name);

	/// Return count of nodes
	inline unsigned int n_nodes() const {
		return n_nodes_;
	}

	/// Return count of elements
	inline unsigned int n_elements() const {
		return n_elements_;
	}

    /**
     *  Reads ElementData sections of opened VTK file.
     *
     *  Implements @p BaseMeshReader::get_element_data.
     */
    template<typename T>
    typename ElementDataCache<T>::ComponentDataPtr get_element_data( std::string field_name, double time, unsigned int n_entities,
    		unsigned int n_components, bool &actual, std::vector<int> const & el_ids, unsigned int component_idx);

protected:
	/// Empty constructor only for tests.
	VtkMeshReader() {}

    /// Reads table of DataArray headers through pugixml interface
    void make_header_table() override;

    /// Helper method that create DataArray header of given xml node (used from \p make_header_table)
    DataArrayAttributes create_header(pugi::xml_node node);

    /// Get DataType by value of string
	DataType get_data_type(std::string type_str);

	/// Return size of value of data_type.
	unsigned int type_value_size(DataType data_type);

	/// Parse ascii data to data cache and return its.
	template<typename T>
	typename ElementDataCache<T>::CacheData parse_ascii_data(unsigned int size_of_cache, unsigned int n_components,
			unsigned int n_entities, std::string data_str);

	/// Parse binary data to data cache and return its.
	template<typename T>
	typename ElementDataCache<T>::CacheData parse_binary_data(unsigned int size_of_cache, unsigned int n_components,
			unsigned int n_entities, unsigned int data_pos, VtkMeshReader::DataType value_type);

	/// Uncompress and parse binary compressed data to data cache and return its.
	template<typename T>
	typename ElementDataCache<T>::CacheData parse_compressed_data(unsigned int size_of_cache, unsigned int n_components,
			unsigned int n_entities, unsigned int data_pos, VtkMeshReader::DataType value_type);

	/// Set count of nodes and elements.
	void read_base_vtk_attributes();

	/// Set @p appended_stream_ for reading AppendedData and find its position in input file
	void set_appended_stream(const FilePath &file_name);

	pugi::xml_document doc_;
    pugi::xml_parse_result parse_result_;

    /// count of nodes
    unsigned int n_nodes_;

    /// count of elements
    unsigned int n_elements_;

    /// header type of VTK file (only for appended data)
    DataType header_type_;

    /// variants of data format (ascii, appended, compressed appended)
    DataFormat data_format_;

    /// Table with data of DataArray headers
    HeaderTable header_table_;

    /// File name (for better error messages)
    std::string f_name_;

    /// input stream allow read appended data, used only if this tag exists
    std::istream *appended_stream_;

    /// position of appended data in file, used only if this tag exists
    unsigned int appended_pos_;

    /// store count of read entities
    unsigned int n_read_;

};

#endif	/* MSH_VTK_READER_HH */

