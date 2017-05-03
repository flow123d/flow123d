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
 * @file    msh_basereader.hh
 * @brief
 * @author  dalibor
 */

#ifndef MSH_BASE_READER_HH
#define	MSH_BASE_READER_HH


#include <string>
#include <vector>
#include <istream>

#include "mesh/element_data_cache.hh"
#include "mesh/mesh.h"
#include "system/system.hh"
#include "system/tokenizer.hh"



/// Types of VTK data (value 'undefined' for empty value)
enum DataType {
	int8, uint8, int16, uint16, int32, uint32, int64, uint64, float32, float64, undefined
};


/***********************************
 * Structure to store the information from a header of \\$ElementData (GMSH file) or DataArray (VTK file) section.
 *
 * Format of GMSH ASCII data sections
 *
   number-of-string-tags (== 2)
     field_name
     interpolation_scheme_name
   number-of-real-tags (==1)
     time_of_dataset
   number-of-integer-tags
     time_step_index (starting from zero)
     number_of_field_components (1, 3, or 9 - i.e. 3d scalar, vector or tensor data)
     number_of entities (nodes or elements)
     partition_index (0 == no partition, not clear if GMSH support reading different partition from different files)
   elm-number value ...
*
*/
struct MeshDataHeader {
    /// Name of field
	std::string field_name;
    /// Currently d'ont used
    std::string interpolation_scheme;
    /// Time of field data (used only for GMSH reader)
    double time;
    /// Currently d'ont used
    unsigned int time_index;
    /// Number of values on one row
    unsigned int n_components;
    /// Number of rows
    unsigned int n_entities;
    /// Currently d'ont used
    unsigned int partition_index;
    /// Position of data in mesh file
    Tokenizer::Position position;
    /// Type of data (used only for VTK reader)
    DataType type;
};


/**
 * Abstract parent of mesh readers.
 *
 * Supported are:
 *  - GMSH reader (class GmshMeshReader)
 *  - VTK reader (class VtkMeshReader)
 */
class BaseMeshReader {
public:
	/// Constructor
	BaseMeshReader(const FilePath &file_name)
	: tok_(file_name) {
		current_cache_ = new ElementDataCacheBase();
	}

	/// Constructor
	BaseMeshReader(std::istream &in)
	: tok_(in) {}

    /**
	 * Reads element data of opened mesh file specified by parameters.
	 *
     * @param field_name field name
     * @param time searched time
     * @param n_entities count of entities (elements)
     * @param n_components count of components (size of returned data is given by n_entities*n_components)
     * @param actual Set to rue if the stream position is just after the header,
     *               set to false either before first header is found or at EOF
     * @param el_ids vector of ids of elements
     * @param component_idx component index of MultiField
	 */
    template<typename T>
    typename ElementDataCache<T>::ComponentDataPtr get_element_data( std::string field_name, double time, unsigned int n_entities,
    		unsigned int n_components, bool &actual, std::vector<int> const & el_ids, unsigned int component_idx)
    { ASSERT(false).error("Method get_element_data must be implement in descendant."); }

protected:
    /**
	 * Find data header for given time and field.
	 */
	virtual MeshDataHeader & find_header(double time, std::string field_name)=0;

    /**
     * Reads table of data headers specific for each descendants.
     */
    virtual void make_header_table()=0;

    /// Cache with last read element data
    ElementDataCacheBase *current_cache_;

    /// Tokenizer used for reading ASCII file format.
    Tokenizer tok_;
};

#endif	/* MSH_BASE_READER_HH */
