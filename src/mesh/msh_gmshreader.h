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
#include "system/tokenizer.hh"

class Mesh;
class GMSH_DataHeader;
class FilePath;

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
    void read_mesh(Mesh* mesh);

    /**
     *  Reads ElementData section of opened GMSH file. Only scalar data are currently supported.
     *  Given @p mesh parameter is used to identify element IDs. The output vector @p data, is
     *  cleared and resized to the size equal to the number of elements in the mesh. Data are stored in the
     *  natural ordering of the ilements in the @p mesh.
     */
    void read_element_data(const std::string &field_name, std::vector<double> &data, Mesh *mesh);

private:
    /**
     * Read section '$PhysicalNames' of the GMSH file and save the physical sections as regions in the RegionDB.
     */
    void read_physical_names(Tokenizer &in);

    /**
     * private method for reading of nodes
     */
    void read_nodes(Tokenizer &in, Mesh*);
    /**
     * private method for reading of elements - in process of implementation
     */
    void read_elements(Tokenizer &in, Mesh*);
    /**
     *
     */
    void read_data_header(Tokenizer &tok, GMSH_DataHeader &head);


    /// Tokenizer used for reading ASCII GMSH file format.
    Tokenizer tok_;
};

#endif	/* _GMSHMESHREADER_H */

