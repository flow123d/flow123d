/*!
 *
 * Copyright (C) 2007 Technical University of Liberec.  All rights reserved.
 *
 * Please make a following refer to Flow123d on your project site if you use the program for any purpose,
 * especially for academic research:
 * Flow123d, Research Centre: Advanced Remedial Technologies, Technical University of Liberec, Czech Republic
 *
 * This program is free software; you can redistribute it and/or modify it under the terms
 * of the GNU General Public License version 3 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 021110-1307, USA.
 *
 *
 * $Id$
 * $Revision$
 * $LastChangedBy$
 * $LastChangedDate$
 *
 *
 * @file   msh_gmshreader.h
 * @author dalibor
 *
 * @date October 3, 2010, 11:23 AM
 */

#ifndef _GMSHMESHREADER_H
#define	_GMSHMESHREADER_H

#include <string>
#include <istream>
#include <vector>
#include <map>


#include "system/tokenizer.hh"
#include "mesh/region.hh"

class Mesh;
class GMSH_DataHeader;
class FilePath;



/***********************************
 * Structure to store the information from a header of $ElementData section.
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
    /// True if the stream position is just after the header.
    /// False either before first header is found or at EOF.
    bool actual;
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
};


class GmshMeshReader {
public:
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
     */
    void read_mesh(Mesh* mesh, const RegionDB::MapElementIDToRegionID *el_to_reg_map=NULL);

    /**
     *  Reads ElementData sections of opened GMSH file. The file is serached for the $ElementData section with header
     *  that match the given @p search_header (same field_name, time of the next section is the first greater then
     *  that given in the @p search_header). If such section has not been yet read, we read the data section into
     *  raw buffer @p data. The map @p id_to_idx is used to convert IDs that marks individual input rows/entities into
     *  indexes to the raw buffer. The buffer must have size at least @p search_header.n_components * @p search_header.n_entities.
     *  Indexes in the map must be smaller then @p search_header.n_entities.
     *
     *  Possible optimizations:
     *  If the map ID lookup seem slow, we may assume that IDs are in increasing order, use simple array of IDs instead of map
     *  and just check that they comes in in correct order.
     */
    void read_element_data( const GMSH_DataHeader &search_header,
            double *data, std::vector<int> const & el_ids);

private:
    /**
     * Read section '$PhysicalNames' of the GMSH file and save the physical sections as regions in the RegionDB.
     *
     * Region Labels starting with '!' are treated as boundary regions. Elements of these regions are used just to
     * assign regions to the boundary and are not used in actual FEM computations.
     */
    void read_physical_names(Tokenizer &in, Mesh * mesh);

    /**
     * private method for reading of nodes
     */
    void read_nodes(Tokenizer &in, Mesh*);
    /**
     * private method for reading of elements - in process of implementation
     */
    void read_elements(Tokenizer &in, Mesh*, const RegionDB::MapElementIDToRegionID *el_to_reg_map=NULL);
    /**
     *
     */
    void read_data_header(Tokenizer &tok, GMSH_DataHeader &head);


    /// Tokenizer used for reading ASCII GMSH file format.
    Tokenizer tok_;
    /// Last read header of ElementData section.
    GMSH_DataHeader last_header;
};

#endif	/* _GMSHMESHREADER_H */

