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

#include "mesh/element_data_cache.hh"
#include "mesh/mesh.h"
#include "system/system.hh"


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
	BaseMeshReader() {}

    /**
     * Reads @p mesh from input file.
     */
    virtual void read_mesh(Mesh* mesh) =0;

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
    /// private method for reading of nodes
    virtual void read_nodes(Mesh*) =0;

    /// Method for reading of elements.
    virtual void read_elements(Mesh*) =0;

    /// Cache with last read element data
    ElementDataCacheBase *current_cache_;

};

#endif	/* MSH_BASE_READER_HH */
