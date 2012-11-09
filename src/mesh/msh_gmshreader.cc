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
 * @file
 * @ingroup mesh
 * @brief
 * @author dalibor
 * 
 * @date Created on October 3, 2010, 11:32 AM
 */

#include "msh_gmshreader.h"
#include "mesh/nodes.hh"
#include <fstream>
#include <string>

#include <boost/tokenizer.hpp>
#include "boost/lexical_cast.hpp"


using namespace std;

Tokenizer::Tokenizer( istream &in)
: in_(in), line_counter_(0), line_tokenizer_(line_,  boost::char_separator<char>("\t \n"))
{}



void Tokenizer::next_line() {
    line_="";
    while ( line_ == "") { std::getline( in_, line_); boost::trim( line_ ); line_counter_++; }
    line_tokenizer_.assign(line_);
    tok_ = line_tokenizer_.begin();
    position = 0;
}

GmshMeshReader::GmshMeshReader()
{
    xprintf(Msg, " - GmshMeshReader()\n");
}



GmshMeshReader::~GmshMeshReader() {
}



void GmshMeshReader::read(const FilePath &file_name, Mesh* mesh) {
    mesh_file = file_name;

    std::ifstream ifs;
    ifs.exceptions ( ifstream::failbit | ifstream::badbit );
    try {
        ifs.open( mesh_file.c_str(), std::ifstream::in );
        read(ifs, mesh);
    }
    catch (ifstream::failure &e)
    {
        /*
         * This doesn't work well, bad bit is set also for failure in getline.
         * 1) check Boost iostreams if thay provides better resolution of exceptions
         * 2) make an open_stream method in FilePath, that can at least check correct oppening
         *    of the stream and turn on the exceptions.
         * 3) If not provided by boost create correct exception types in FilePath
         *    for reporting type of io problem and possibly information about filename and
         *    line.
         */
        if (ifs.bad())
            xprintf(UsrErr,"Can not open GMSH input file: %s\n", mesh_file.c_str());
        if (ifs.fail())
            xprintf(UsrErr,"Can not read GMSH input file: %s. Should be text file.\n", mesh_file.c_str());
    }

    mesh_file ="";
}




void GmshMeshReader::read(istream &in, Mesh *mesh) {
    F_ENTRY;

    ASSERT( mesh , "Argument mesh is NULL.\n");
    INPUT_CHECK( ! in.fail(), "Can not open GMSH input file: %s\n", mesh_file.c_str());
    read_nodes(in, mesh);
    read_elements(in, mesh);

    mesh->setup_topology();
}




void GmshMeshReader::read_nodes(istream &in, Mesh* mesh) {
    using namespace boost;
    xprintf(Msg, " - Reading nodes...");

    skip_to(in, "$Nodes");
    Tokenizer tok(in);
    try {
        tok.next_line();
        DBGMSG("%d\n",tok.line_num());
        unsigned int n_nodes = lexical_cast<unsigned int> (*tok);;

        INPUT_CHECK( n_nodes > 0, "Zero number of nodes in the mesh file %s.", mesh_file.c_str() ); // should throw and catch at level where we know the file name
        mesh->node_vector.reserve(n_nodes);

        for (int i = 0; i < n_nodes; ++i) {
            tok.next_line();

            unsigned int id = lexical_cast<unsigned int> (*tok); ++tok;
            NodeFullIter node = mesh->node_vector.add_item(id);

            node->point()(0)=lexical_cast<double> (*tok); ++tok;
            node->point()(1)=lexical_cast<double> (*tok); ++tok;
            node->point()(2)=lexical_cast<double> (*tok);
        }

    } catch (bad_lexical_cast &) {
        xprintf(UsrErr, "Wrong number at line %d of the '$Nodes' section in mesh file '%s'\n", tok.line_num(), mesh_file.c_str());
    }
    xprintf(Msg, " %d nodes read. \n ", mesh->node_vector.size());
}





void GmshMeshReader::read_elements(istream &in, Mesh * mesh) {
    using namespace boost;
    xprintf(Msg, " - Reading elements...");

    skip_to(in, "$Elements");
    Tokenizer tok(in);

    try {
        tok.next_line();
        unsigned int n_elements = lexical_cast<unsigned int> (*tok);
        INPUT_CHECK( n_elements > 0, "Zero number of elements in the mesh file %s.\n", mesh_file.c_str());
        mesh->element.reserve(n_elements);

        for (unsigned int i = 0; i < n_elements; ++i) {
            tok.next_line();

            unsigned int id = lexical_cast<unsigned int>(*tok); ++tok;

            ElementFullIter ele(mesh->element.add_item(id));

            //get element type: supported:
            //  1 Line (2 nodes)
            //  2 Triangle (3 nodes)
            //  4 Tetrahedron (4 nodes)
            unsigned int type = lexical_cast<unsigned int>(*tok); ++tok;
            switch (type) {
                case 1:
                    ele->dim_ = 1;
                    break;
                case 2:
                    ele->dim_ = 2;
                    break;
                case 4:
                    ele->dim_ = 3;
                    break;
                default:
                    xprintf(UsrErr, "Element %d is of the unsupported type %d\n", ele.id(), type);
            }

            //get number of tags (at least 2)
            unsigned int n_tags = lexical_cast<unsigned int>(*tok); ++tok;
            INPUT_CHECK(n_tags >= 2, "At least two element tags have to be defined for elment with id=%d\n", id);
            //get tags 1 and 2
            ele->mid = lexical_cast<unsigned int>(*tok); ++tok;
            unsigned int rid = lexical_cast<unsigned int>(*tok); ++tok; // region number, we do not store this
            //get remaining tags
            if (n_tags > 2)  { ele->pid = lexical_cast<unsigned int>(*tok); ++tok; } // save partition number in the new GMSH format
            for (unsigned int ti = 3; ti < n_tags; ti++) ++tok;         //skip remaining tags

            // allocate element arrays TODO: should be in mesh class
            ele->node = new Node * [ele->n_nodes()];
            ele->edges_ = new Edge * [ele->n_sides()];
            ele->boundaries_ = new Boundary * [ele->n_sides()];
            ele->mesh_ = mesh;

            FOR_ELEMENT_SIDES(ele, si) {
                ele->edges_[ si ]=NULL;
                ele->boundaries_[si] =NULL;
            }

            unsigned int ni;
            FOR_ELEMENT_NODES(ele, ni) {
                unsigned int node_id = lexical_cast<unsigned int>(*tok); ++tok;
                NodeIter node = mesh->node_vector.find_id( node_id );
                INPUT_CHECK( node!=mesh->node_vector.end() , "Unknown node id %d in specification of element with id=%d.\n", node_id, id);
                ele->node[ni] = node;
            }
        }

    } catch (bad_lexical_cast &) {
        xprintf(UsrErr, "Wrong number at line %d of the '$Elements' section in mesh file '%s'\n", tok.line_num(), mesh_file.c_str());
    }

    xprintf(Msg, " %d elements read. \n", mesh->n_elements());
}
