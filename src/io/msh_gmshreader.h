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

#include <boost/exception/info.hpp>  // for error_info::~error_info<Tag, T>
#include <map>                       // for map, map<>::value_compare
#include <string>                    // for string
#include <vector>                    // for vector
#include "io/msh_basereader.hh"      // for MeshDataHeader, BaseMeshReader
#include "system/exceptions.hh"      // for ExcStream, operator<<, EI, TYPED...

class ElementDataCacheBase;
class FilePath;
class Mesh;



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
     * Construct the GMSH format reader from given FilePath.
     * This opens the file for reading.
     */
    GmshMeshReader(const FilePath &file_name);

    /**
     * Destructor close the file if opened.
     */
    virtual ~GmshMeshReader();

    /**
     * Read section '$PhysicalNames' of the GMSH file and save the physical sections as regions in the RegionDB.
     *
     * Region Labels starting with '!' are treated as boundary regions. Elements of these regions are used just to
     * assign regions to the boundary and are not used in actual FEM computations.
     */
    void read_physical_names(Mesh * mesh) override;

    /**
     * Empty method for GMSH reader now.
     *
     * Implements @p BaseMeshReader::check_compatible_mesh.
     *
     * OBSOLETE method - will be replace with Mesh::check_compatible_mesh after merge fields!
     */
    void check_compatible_mesh(Mesh &mesh) override;

    /**
     * Finds GMSH data header for ElementData given by time and field_name and stores it as \p actual_header_.
     */
    MeshDataHeader & find_header(HeaderQuery &header_query) override;

protected:
	/**
	 * Map of ElementData sections in GMSH file.
	 *
	 * For each field_name contains vector of MeshDataHeader.
	 * Headers are sorted by time in ascending order.
	 */
	typedef typename std::map< std::string, std::vector<MeshDataHeader> > HeaderTable;

    /**
     * private method for reading of nodes
     */
    void read_nodes(Mesh * mesh);
    /**
     *  Method for reading of elements.
     *  Input of the mesh allows changing regions within the input CON file.
     * Read section '$PhysicalNames' of the GMSH file and save the physical sections to general data structure.
     *
     * Region Labels starting with '!' are treated as boundary regions. Elements of these regions are used just to
     * assign regions to the boundary and are not used in actual FEM computations.
     */
    void read_elements(Mesh * mesh);
    /**
     * Reads the header from the tokenizer @p tok and return it as the second parameter.
     */
    void read_data_header(MeshDataHeader &head);
    /**
     * Reads table of ElementData headers from the tokenizer file.
     */
    void make_header_table() override;
    /**
     * Implements @p BaseMeshReader::read_element_data.
     */
    void read_element_data(ElementDataCacheBase &data_cache, MeshDataHeader actual_header, unsigned int n_components,
    		bool boundary_domain) override;


    /// Table with data of ElementData headers
    HeaderTable header_table_;
};

#endif	/* _GMSHMESHREADER_H */


