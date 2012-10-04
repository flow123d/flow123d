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

#include "mesh/mesh.h"
#include <string>
#include <istream>
#include <boost/tokenizer.hpp>


class Tokenizer {
public:
    typedef boost::tokenizer<boost::char_separator<char> > BT;

    Tokenizer( istream &in);
    void next_line();
    inline const std::string & operator *() const
        { return *tok_; }
    inline BT::iterator & operator ++()
        {position++; return ++tok_;}
    inline unsigned int line_num() const
        {return line_counter_;}
private:
    istream &in_;
    string line_;
    unsigned int line_counter_;

    unsigned int position;
    BT::iterator tok_;
    BT line_tokenizer_;
};



class GmshMeshReader {
private:
    std::string mesh_file;

    /**
     * private method for reading of nodes
     */
    void read_nodes(istream &in, Mesh*);
    /**
     * private method for reading of elements - in process of implementation
     */
    void read_elements(istream &in, Mesh*);

public:
    GmshMeshReader();
    ~GmshMeshReader();

    /**
     *  Reads @p mesh from file with given @p file_name.
     */
    void read(const FilePath &file_name, Mesh* mesh);

    /**
     *  Reads @p mesh from given stream @p in.
     */
    void read(istream &in, Mesh *mesh);

};

#endif	/* _GMSHMESHREADER_H */

