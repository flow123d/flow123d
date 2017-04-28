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
 * @file    msh_gmshreader.h
 * @brief   
 * @author  dalibor
 */

#ifndef _GMSHMESHREADER_H
#define	_GMSHMESHREADER_H

#include <string>
#include <istream>
#include <vector>
#include <map>


#include "system/tokenizer.hh"
#include "mesh/region.hh"
#include "mesh/element_data_cache.hh"
#include "mesh/msh_basereader.hh"
#include "input/input_exception.hh"

class Mesh;
class FilePath;



/***********************************
 * Structure to store the information from a header of \\$ElementData section.
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

struct GMSH_DataHeader {
    std::string field_name;
    /// Currently ont used
    std::string interpolation_scheme;
    double time;
    /// Currently ont used
    unsigned int time_index;
    /// Number of values on one row
    unsigned int n_components;
    /// Number of rows
    unsigned int n_entities;
    /// ?? Currently ont used
    unsigned int partition_index;
    /// Position of data in mesh file
    Tokenizer::Position position;
};


class GmshMeshReader : public BaseMeshReader {
public:
	TYPEDEF_ERR_INFO(EI_FieldName, std::string);
	TYPEDEF_ERR_INFO(EI_GMSHFile, std::string);
	TYPEDEF_ERR_INFO(EI_Time, double);
	TYPEDEF_ERR_INFO(EI_Type, std::string);
	TYPEDEF_ERR_INFO(EI_TokenizerMsg, std::string);
	TYPEDEF_ERR_INFO(EI_Section, std::string);
	TYPEDEF_ERR_INFO(EI_ElementId, int);
	TYPEDEF_ERR_INFO(EI_ElementType, int);
	DECLARE_INPUT_EXCEPTION(ExcFieldNameNotFound,
			<< "No data for field: "<< EI_FieldName::qval
			<< " and time: "<< EI_Time::val
			<< " in the input file: "<< EI_GMSHFile::qval);
	DECLARE_EXCEPTION(ExcWrongFormat,
			<< "Wrong format of " << EI_Type::val << ", " << EI_TokenizerMsg::val << "\n"
			<< "in the input file: " << EI_GMSHFile::qval);
	DECLARE_EXCEPTION(ExcMissingSection,
			<< "Missing section " << EI_Section::qval << " in the GMSH input file: " << EI_GMSHFile::qval);
	DECLARE_EXCEPTION(ExcUnsupportedType,
			<< "Element " << EI_ElementId::val << "in the GMSH input file " << EI_GMSHFile::qval
			<< " is of the unsupported type " << EI_ElementType::val );

	/**
	 * Map of ElementData sections in GMSH file.
	 *
	 * For each field_name contains vector of GMSH_DataHeader.
	 * Headers are sorted by time in ascending order.
	 */
	typedef typename std::map< std::string, std::vector<GMSH_DataHeader> > HeaderTable;

    /**
     * Construct the GMSH format reader from given filename.
     * This opens the file for reading.
     */
    GmshMeshReader(const FilePath &file_name);
    /**
     * Construct the GMSH format reader from given input stream.
     * The input stream should be correctly opened. To get correct information about
     * line numbers there should be no previous reading from the stream.
     */
    GmshMeshReader(std::istream &in);

    /**
     * Destructor close the file if opened.
     */
    ~GmshMeshReader();

    /**
     *  Reads @p mesh from the GMSH file.
     *  Input of the mesh allows changing regions within the input CON file.
     *
     *  Implements @p BaseMeshReader::read_mesh.
     */
    void read_mesh(Mesh* mesh);

    /**
     * Read section '$PhysicalNames' of the GMSH file and save the physical sections as regions in the RegionDB.
     *
     * Region Labels starting with '!' are treated as boundary regions. Elements of these regions are used just to
     * assign regions to the boundary and are not used in actual FEM computations.
     */
    void read_physical_names(Mesh * mesh);

    /**
     *  Reads ElementData sections of opened GMSH file. The file is searched for the \\$ElementData section with header
     *  that match the given @p field_name and @p time of the next section is the first greater then that given in input
     *  parameters). If such section has not been yet read, we read the data section into raw buffer @p data. The map
     *  @p id_to_idx is used to convert IDs that marks individual input rows/entities into indexes to the raw buffer.
     *  The buffer must have size at least @p n_components * @p n_entities. Indexes in the map must be smaller then
     *  @p n_entities.
     *  If the @p data buffer is updated we set actual to true.
     *
     *  Possible optimizations:
     *  If the map ID lookup seem slow, we may assume that IDs are in increasing order, use simple array of IDs instead of map
     *  and just check that they comes in in correct order.
     *
     *  Implements @p BaseMeshReader::get_element_data.
     */
    template<typename T>
    typename ElementDataCache<T>::ComponentDataPtr get_element_data( std::string field_name, double time, unsigned int n_entities,
    		unsigned int n_components, bool &actual, std::vector<int> const & el_ids, unsigned int component_idx);

protected:
    /**
     * private method for reading of nodes
     */
    void read_nodes(Mesh*);
    /**
     *  Method for reading of elements.
     *  Input of the mesh allows changing regions within the input CON file.
     */
    void read_elements(Mesh*);
    /**
     * Reads the header from the tokenizer @p tok and return it as the second parameter.
     */
    void read_data_header(GMSH_DataHeader &head);
    /**
     * Reads table of ElementData headers from the tokenizer file.
     */
    void make_header_table() override;
    /**
     * Finds GMSH data header for ElementData given by time and field_name and return it as the first parameter.
     */
    GMSH_DataHeader & find_header(double time, std::string field_name);


    /// Tokenizer used for reading ASCII GMSH file format.
    Tokenizer tok_;
    /// Table with data of ElementData headers
    HeaderTable header_table_;
};

#endif	/* _GMSHMESHREADER_H */


