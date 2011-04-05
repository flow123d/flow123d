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
 * @file   gmshmeshreader.h
 * @author dalibor
 *
 * @date October 3, 2010, 11:23 AM
 */

#ifndef _GMSHMESHREADER_H
#define	_GMSHMESHREADER_H

#include "msh_reader.h"
#include "mesh.h"

class GmshMeshReader : public MeshReader {
private:
    void read_nodes(FILE*, Mesh*);
    void read_elements(FILE*, Mesh*);

    char supported_element_type(int);
    void parse_element_line(ElementVector&, char*, Mesh* mesh);
    void element_type_specific(ElementFullIter);
    void element_allocation_independent(ElementFullIter);

public:
    GmshMeshReader();
    ~GmshMeshReader();

    /**
     *  Read mesh from file
     */
    void read(const std::string &file_name, Mesh* mesh);

};

#endif	/* _GMSHMESHREADER_H */

