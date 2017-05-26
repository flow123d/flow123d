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
#include <vector>
#include <map>


#include "mesh/region.hh"
#include "io/element_data_cache.hh"
#include "io/msh_basereader.hh"
#include "input/input_exception.hh"

class Mesh;
class FilePath;



class GmshMeshReader : public BaseMeshReader {
public:
	TYPEDEF_ERR_INFO(EI_GMSHFile, std::string);
	TYPEDEF_ERR_INFO(EI_Section, std::string);
	TYPEDEF_ERR_INFO(EI_ElementId, int);
	TYPEDEF_ERR_INFO(EI_ElementType, int);
	DECLARE_EXCEPTION(ExcMissingSection,
			<< "Missing section " << EI_Section::qval << " in the GMSH input file: " << EI_GMSHFile::qval);
	DECLARE_EXCEPTION(ExcUnsupportedType,
			<< "Element " << EI_ElementId::val << "in the GMSH input file " << EI_GMSHFile::qval
			<< " is of the unsupported type " << EI_ElementType::val );

	/**
	 * Map of ElementData sections in GMSH file.
	 *
	 * For each field_name contains vector of MeshDataHeader.
	 * Headers are sorted by time in ascending order.
	 */
	typedef typename std::map< std::string, std::vector<MeshDataHeader> > HeaderTable;

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
     * Empty method for GMSH reader now.
     *
     * Implements @p BaseMeshReader::check_compatible_mesh.
     */
    void check_compatible_mesh(Mesh &mesh) override;

    /**
     * Read section '$Nodes' of the GMSH file and save the physical sections to general data structure.
     *
     * Implements @p BaseMeshReader::read_nodes_data.
     */
    NodeDataTable read_nodes_data() override;

    /**
     * Read section '$Elements' of the GMSH file and save the physical sections to general data structure.
     */
    ElementDataTable read_elements_data();

    /**
     * Read section '$PhysicalNames' of the GMSH file and save the physical sections to general data structure.
     *
     * Region Labels starting with '!' are treated as boundary regions. Elements of these regions are used just to
     * assign regions to the boundary and are not used in actual FEM computations.
     */
    PhysicalNamesDataTable read_physical_names_data();

protected:
    /**
     * Reads the header from the tokenizer @p tok and return it as the second parameter.
     */
    void read_data_header(MeshDataHeader &head);
    /**
     * Reads table of ElementData headers from the tokenizer file.
     */
    void make_header_table() override;
    /**
     * Finds GMSH data header for ElementData given by time and field_name and return it as the first parameter.
     */
    MeshDataHeader & find_header(double time, std::string field_name) override;
    /**
     * Implements @p BaseMeshReader::read_element_data.
     */
    void read_element_data(ElementDataCacheBase &data_cache, MeshDataHeader actual_header, unsigned int size_of_cache,
    		unsigned int n_components, std::vector<int> const & el_ids) override;

    /// Implements @p BaseMeshReader::data_section_name.
    std::string data_section_name() override {
    	return "$ElementData";
    }


    /// Table with data of ElementData headers
    HeaderTable header_table_;
};

#endif	/* _GMSHMESHREADER_H */


