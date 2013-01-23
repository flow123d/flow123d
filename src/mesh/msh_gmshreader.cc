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

#include "mesh/mesh.h"
#include "mesh/nodes.hh"
#include <istream>
#include <string>

#include "system/system.hh"
#include "system/tokenizer.hh"
#include "boost/lexical_cast.hpp"


using namespace std;


GmshMeshReader::GmshMeshReader(const FilePath &file_name)
: tok_(file_name)
{}



GmshMeshReader::GmshMeshReader(std::istream &in)
: tok_(in)
{}



GmshMeshReader::~GmshMeshReader()   // Tokenizer close the file automatically
{}



void GmshMeshReader::read_mesh(Mesh* mesh) {
    F_ENTRY;

    ASSERT( mesh , "Argument mesh is NULL.\n");
    read_physical_names(tok_);
    read_nodes(tok_, mesh);
    read_elements(tok_, mesh);
    mesh->setup_topology();
}



void GmshMeshReader::read_nodes(Tokenizer &tok, Mesh* mesh) {
    using namespace boost;
    xprintf(Msg, "- Reading nodes...");

    if (! tok.skip_to("$Nodes")) xprintf(UsrErr,"Missing section '$Nodes' in the GMSH input file: %s\n",tok.f_name().c_str());
    try {
        tok.next_line(false);
        unsigned int n_nodes = lexical_cast<unsigned int> (*tok);;
        INPUT_CHECK( n_nodes > 0, "Zero number of nodes, %s.\n", tok.position_msg().c_str() );
        ++tok; // end of line

        mesh->node_vector.reserve(n_nodes);
        for (int i = 0; i < n_nodes; ++i) {
            tok.next_line();

            unsigned int id = lexical_cast<unsigned int> (*tok); ++tok;
            NodeFullIter node = mesh->node_vector.add_item(id);

            node->point()(0)=lexical_cast<double> (*tok); ++tok;
            node->point()(1)=lexical_cast<double> (*tok); ++tok;
            node->point()(2)=lexical_cast<double> (*tok); ++tok;
            ++tok; // skip mesh size parameter
        }

    } catch (bad_lexical_cast &) {
        xprintf(UsrErr, "Wrong format of number, %s.\n", tok.position_msg().c_str());
    }
    xprintf(Msg, " %d nodes read. \n", mesh->node_vector.size());
}



void GmshMeshReader::read_elements(Tokenizer &tok, Mesh * mesh) {
    using namespace boost;
    xprintf(Msg, "- Reading elements...");

    if (! tok.skip_to("$Elements")) xprintf(UsrErr,"Missing section '$Elements' in the GMSH input file: %s\n",tok.f_name().c_str());
    try {
        tok.next_line(false);
        unsigned int n_elements = lexical_cast<unsigned int> (*tok);
        INPUT_CHECK( n_elements > 0, "Zero number of elements, %s.\n", tok.position_msg().c_str());
        ++tok; // end of line

        mesh->element.reserve(n_elements);

        for (unsigned int i = 0; i < n_elements; ++i) {
            tok.next_line();

            unsigned int id = lexical_cast<unsigned int>(*tok); ++tok;


            //get element type: supported:
            //  1 Line (2 nodes)
            //  2 Triangle (3 nodes)
            //  4 Tetrahedron (4 nodes)
            unsigned int type = lexical_cast<unsigned int>(*tok); ++tok;
            unsigned int dim;
            switch (type) {
                case 1:
                    dim = 1;
                    break;
                case 2:
                    dim = 2;
                    break;
                case 4:
                    dim = 3;
                    break;
                default:
                    xprintf(UsrErr, "Element %d is of the unsupported type %d\n", id, type);
            }

            //get number of tags (at least 2)
            unsigned int n_tags = lexical_cast<unsigned int>(*tok);
            INPUT_CHECK(n_tags >= 2, "At least two element tags have to be defined for element with id=%d, %s.\n",
                    id, tok.position_msg().c_str());
            ++tok;

            //get tags 1 and 2
            unsigned int region_id = lexical_cast<unsigned int>(*tok); ++tok;
            unsigned int object_id = lexical_cast<unsigned int>(*tok); ++tok; // GMSH region number, we do not store this
            //get remaining tags
            unsigned int partition_id;
            if (n_tags > 2)  { partition_id = lexical_cast<unsigned int>(*tok); ++tok; } // save partition number from the new GMSH format
            for (unsigned int ti = 3; ti < n_tags; ti++) ++tok;         //skip remaining tags

            // allocate element arrays TODO: should be in mesh class
            Element *ele;
            Region region_idx = Region::db().add_region( region_id, dim );
            if (region_idx.is_boundary()) {
                ele = mesh->bc_elements.add_item(id);
            } else {
                ele = mesh->element.add_item(id);
            }
            ele->dim_=dim;
            ele->region_=region_idx;
            ele->pid=partition_id;

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
                unsigned int node_id = lexical_cast<unsigned int>(*tok);
                NodeIter node = mesh->node_vector.find_id( node_id );
                INPUT_CHECK( node!=mesh->node_vector.end() ,
                        "Unknown node id %d in specification of element with id=%d, %s.\n",
                        node_id, id, tok.position_msg().c_str());
                ele->node[ni] = node;
                ++tok;
            }
        }

    } catch (bad_lexical_cast &) {
        xprintf(UsrErr, "Wrong format of number, %s.\n", tok.position_msg().c_str());
    }

    xprintf(Msg, " %d elements read. \n", mesh->n_elements());
}



void GmshMeshReader::read_physical_names(Tokenizer &tok) {
    using namespace boost;

    if (! tok.skip_to("$PhysicalNames", "$Nodes") ) return;
    try {
        tok.next_line(false);
        unsigned int n_physicals = lexical_cast<unsigned int> (*tok);
        ++tok; // end of line

        for (unsigned int i = 0; i < n_physicals; ++i) {
            tok.next_line();
            // format of one line:
            // dim    physical-id    physical-name

            unsigned int dim = lexical_cast<unsigned int>(*tok); ++tok;
            unsigned int id = lexical_cast<unsigned int>(*tok); ++tok;
            string name = *tok; ++tok;

            bool boundary =  ( name.size() != 0 && name[0] == '.' );
            Region::db().add_region(id, name, dim, boundary);
        }

    } catch (bad_lexical_cast &) {
        xprintf(UsrErr, "Wrong format of number, %s.\n", tok.position_msg().c_str());
    }
}


/***********************************
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

*/


/// Structure to which store the information form the header of data section.
///
struct GMSH_DataHeader {
    string field_name;
    string interpolation_scheme;
    double time;
    unsigned int time_index;
    unsigned int n_components;
    unsigned int n_entities;
    unsigned int partition_index;
};

// Is assumed to be called just after tok.skip_to("..")
// reads the header from the tokenizer @p tok and return it as the second parameter
void GmshMeshReader::read_data_header(Tokenizer &tok, GMSH_DataHeader &head) {
    using namespace boost;

    // string tags
    tok.next_line();
    unsigned int n_str = lexical_cast<unsigned int>(*tok);
    head.field_name="";
    head.interpolation_scheme = "";
    if (n_str > 0) {
        tok.next_line(); n_str--;
        head.field_name= *tok; //  unquoted by tokenizer if needed
    }
    if (n_str > 0) {
        tok.next_line(); n_str--;
        head.interpolation_scheme = *tok;
    }
    for(;n_str>0;n_str--) tok.next_line(); // skip possible remaining tags

    //real tags
    tok.next_line();
    unsigned int n_real = lexical_cast<unsigned int>(*tok);
    head.time=0.0;
    if (n_real>0) {
        tok.next_line(); n_real--;
        head.time=lexical_cast<double>(*tok);
    }
    for(;n_real>0;n_real--) tok.next_line();

    // int tags
    tok.next_line();
    unsigned int n_int = lexical_cast<unsigned int>(*tok);
    head.time_index=0;
    head.n_components=1;
    head.n_entities=0;
    head.partition_index=0;
    if (n_int>0) {
        tok.next_line(); n_int--;
        head.time_index=lexical_cast<unsigned int>(*tok);
    }
    if (n_int>0) {
        tok.next_line(); n_int--;
        head.n_components=lexical_cast<unsigned int>(*tok);
    }
    if (n_int>0) {
        tok.next_line(); n_int--;
        head.n_entities=lexical_cast<unsigned int>(*tok);
    }
    for(;n_int>0;n_int--) tok.next_line();
}



void GmshMeshReader::read_element_data( const  string &field_name, vector<double> &data, Mesh *mesh) {
    using namespace boost;

    // find the data section
    GMSH_DataHeader head;
    head.field_name="";
    while (! tok_.eof() && head.field_name!=field_name) {
        tok_.skip_to("$ElementData");
        read_data_header(tok_, head);
    }

    INPUT_CHECK(head.field_name==field_name, "Can not find data field '%s' in input file '%s'.\n",
            field_name.c_str(), tok_.f_name().c_str());
    ASSERT( mesh , "Argument mesh is NULL.\n");

    // read the data
    vector<double>(data).swap(data); // clear vector by swapping
    data.reserve(head.n_entities);
    int id;
    double scalar_value;
    for (unsigned int i = 0; i < head.n_entities; ++i) {
        tok_.next_line();
        id = lexical_cast<unsigned int>(*tok_); ++tok_;
        scalar_value = lexical_cast<double>(*tok_); ++tok_;

        ElementFullIter ele = mesh->element.find_id(id);
        data[ ele.index() ] = scalar_value;
    }

    xprintf(Msg, " %d values of pressure read. O.K.\n", data.size());
}

