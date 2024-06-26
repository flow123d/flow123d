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
 * @file    msh_gmshreader.cc
 * @ingroup mesh
 * @brief   
 * @author  dalibor
 */

#include <istream>
#include <string>
#include <limits>

#include "msh_gmshreader.h"
#include "io/element_data_cache_base.hh"

#include "system/system.hh"
#include "system/tokenizer.hh"
#include "boost/lexical_cast.hpp"

#include "mesh/mesh.h"
#include "mesh/bc_mesh.hh"



using namespace std;


GmshMeshReader::GmshMeshReader(const FilePath &file_name)
: BaseMeshReader(file_name)
{
    tok_.set_comment_pattern( "#");
    data_section_name_ = "$ElementData";
    has_compatible_mesh_ = false;
    make_header_table();
}



GmshMeshReader::~GmshMeshReader()   // Tokenizer close the file automatically
{}



void GmshMeshReader::read_nodes(Mesh * mesh) {
    using namespace boost;
    unsigned int n_nodes;
    MessageOut() << "- Reading nodes...";
    tok_.set_position( Tokenizer::Position() );

    if (! tok_.skip_to("$Nodes")) THROW(ExcMissingSection() << EI_Section("$Nodes") << EI_GMSHFile(tok_.f_name()) );
    try {
    	tok_.next_line(false);
        n_nodes = lexical_cast<unsigned int> (*tok_);
        mesh->init_node_vector( n_nodes );
        if (n_nodes == 0) THROW( ExcZeroNodes() << EI_Position(tok_.position_msg()) );
        ++tok_; // end of line

        for (unsigned int i = 0; i < n_nodes; ++i) {
        	tok_.next_line();

            unsigned int id = lexical_cast<unsigned int> (*tok_); ++tok_; // node id
        	arma::vec3 coords;                                         // node coordinates
        	coords(0) = lexical_cast<double> (*tok_); ++tok_;
        	coords(1) = lexical_cast<double> (*tok_); ++tok_;
        	coords(2) = lexical_cast<double> (*tok_); ++tok_;
            ++tok_; // skip mesh size parameter

            mesh->add_node(id, coords);
        }

    } catch (bad_lexical_cast &) {
    	THROW(ExcWrongFormat() << EI_Type("number") << EI_TokenizerMsg(tok_.position_msg()) << EI_MeshFile(tok_.f_name()) );
    }
    MessageOut().fmt("... {} nodes read. \n", n_nodes);
}


void GmshMeshReader::read_elements(Mesh * mesh) {
    using namespace boost;
    MessageOut() << "- Reading elements...";

    if (! tok_.skip_to("$Elements")) THROW(ExcMissingSection() << EI_Section("$Elements") << EI_GMSHFile(tok_.f_name()) );
    try {
    	tok_.next_line(false);
        unsigned int n_elements = lexical_cast<unsigned int> (*tok_);
        if (n_elements == 0) THROW( ExcZeroElements() << EI_Position(tok_.position_msg()) );
        ++tok_; // end of line

        std::vector<unsigned int> node_ids; //node_ids of elements
        node_ids.resize(4); // maximal count of nodes

        mesh->init_element_vector(n_elements);

        for (unsigned int i = 0; i < n_elements; ++i) {
        	tok_.next_line();
            unsigned int id = lexical_cast<unsigned int>(*tok_); ++tok_;

            //get element type: supported:
            //  1 Line (2 nodes)
            //  2 Triangle (3 nodes)
            //  4 Tetrahedron (4 nodes)
            // 15 Point (1 node)
            unsigned int type = lexical_cast<unsigned int>(*tok_); ++tok_;
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
                case 15:
                    dim = 0;
                    break;
                default:
                    dim = 0;
                    THROW(ExcUnsupportedType() << EI_ElementId(id) << EI_ElementType(type) << EI_GMSHFile(tok_.f_name()) );
                    break;
            }

            //get number of tags (at least 2)
            unsigned int n_tags = lexical_cast<unsigned int>(*tok_);
            if (n_tags < 2) THROW( ExcTooManyElementTags() << EI_ElementId(id) << EI_Position(tok_.position_msg()) );
            ++tok_;

            //get tags 1 and 2
            unsigned int region_id = lexical_cast<unsigned int>(*tok_); ++tok_; // region_id
            lexical_cast<unsigned int>(*tok_); ++tok_; // GMSH region number, we do not store this
            //get remaining tags
            unsigned int partition_id = 0;
            if (n_tags > 2)  { partition_id = lexical_cast<unsigned int>(*tok_); ++tok_; } // save partition number from the new GMSH format
            for (unsigned int ti = 3; ti < n_tags; ti++) ++tok_;         //skip remaining tags

            for (unsigned int ni=0; ni<dim+1; ++ni) { // read node ids
            	node_ids[ni] = lexical_cast<unsigned int>(*tok_);
                ++tok_;
            }
            mesh->add_element(id, dim, region_id, partition_id, node_ids);
        }

    } catch (bad_lexical_cast &) {
    	THROW(ExcWrongFormat() << EI_Type("number") << EI_TokenizerMsg(tok_.position_msg()) << EI_MeshFile(tok_.f_name()) );
    }

    MessageOut().fmt("... {} bulk elements, {} boundary elements. \n", mesh->n_elements(), mesh->bc_mesh()->n_elements());
}



void GmshMeshReader::read_physical_names(Mesh * mesh) {
	ASSERT_PTR(mesh).error("Argument mesh is NULL.\n");

    using namespace boost;

    if (! tok_.skip_to("$PhysicalNames", "$Nodes") ) return;
    try {
    	tok_.next_line(false);
        unsigned int n_physicals = lexical_cast<unsigned int> (*tok_);
        ++tok_; // end of line

        for (unsigned int i = 0; i < n_physicals; ++i) {
        	tok_.next_line();
            // format of one line:
            // dim    physical-id    physical-name

            unsigned int dim = lexical_cast<unsigned int>(*tok_); ++tok_;
            unsigned int id = lexical_cast<unsigned int>(*tok_); ++tok_;
            string name = *tok_; ++tok_;
            mesh->add_physical_name( dim, id, name );
        }

    } catch (bad_lexical_cast &) {
    	THROW(ExcWrongFormat() << EI_Type("number") << EI_TokenizerMsg(tok_.position_msg()) << EI_MeshFile(tok_.f_name()) );
    }

}


// Is assumed to be called just after tok.skip_to("..")
// reads the header from the tokenizer @p tok and return it as the second parameter
void GmshMeshReader::read_data_header(MeshDataHeader &head) {
    using namespace boost;
    try {
        // string tags
    	tok_.next_line(false);
        unsigned int n_str = lexical_cast<unsigned int>(*tok_); ++tok_;
        head.field_name="";
        head.interpolation_scheme = "";
        if (n_str > 0) {
        	tok_.next_line(); n_str--;
            head.field_name= *tok_; ++tok_; //  unquoted by tokenizer if needed
        }
        if (n_str > 0) {
        	tok_.next_line(); n_str--;
            head.interpolation_scheme = *tok_; ++tok_;
        }
        for(;n_str>0;n_str--) tok_.next_line(false); // skip possible remaining tags

        //real tags
        tok_.next_line();
        unsigned int n_real = lexical_cast<unsigned int>(*tok_); ++tok_;
        head.time=0.0;
        if (n_real>0) {
        	tok_.next_line(); n_real--;
            head.time=lexical_cast<double>(*tok_); ++tok_;
        }
        for(;n_real>0;n_real--) tok_.next_line(false);

        // int tags
        tok_.next_line();
        unsigned int n_int = lexical_cast<unsigned int>(*tok_); ++tok_;
        head.time_index=0;
        head.n_components=1;
        head.n_entities=0;
        head.partition_index=0;
        if (n_int>0) {
        	tok_.next_line(); n_int--;
            head.time_index=lexical_cast<unsigned int>(*tok_); ++tok_;
        }
        if (n_int>0) {
        	tok_.next_line(); n_int--;
            head.n_components=lexical_cast<unsigned int>(*tok_); ++tok_;
        }
        if (n_int>0) {
        	tok_.next_line(); n_int--;
            head.n_entities=lexical_cast<unsigned int>(*tok_); ++tok_;
        }
        for(;n_int>0;n_int--) tok_.next_line(false);
        head.position = tok_.get_position();
        head.discretization = OutputTime::DiscreteSpace::ELEM_DATA;
    } catch (bad_lexical_cast &) {
    	THROW(ExcWrongFormat() << EI_Type("$ElementData header") << EI_TokenizerMsg(tok_.position_msg()) << EI_MeshFile(tok_.f_name()) );
    }
}



void GmshMeshReader::read_element_data(ElementDataCacheBase &data_cache, MeshDataHeader header) {
    static int imax = std::numeric_limits<int>::max();
    unsigned int id, i_row;
    unsigned int n_bulk_read = 0, n_bdr_read = 0;
    std::vector<int> bulk_el_ids = this->get_element_ids(false); // bulk
    bulk_el_ids.push_back( imax ); // put 'save' item at the end of vector
    vector<int>::const_iterator bulk_id_iter = bulk_el_ids.begin();
    std::vector<int> bdr_el_ids = this->get_element_ids(true); // boundary
    bdr_el_ids.push_back( imax ); // put 'save' item at the end of vector
    vector<int>::const_iterator bdr_id_iter = bdr_el_ids.begin();

    // read @p data buffer as we have correct header with already passed time
    // we assume that @p data buffer is big enough
    tok_.set_position(header.position);

    // read data
    for (i_row = 0; i_row < header.n_entities; ++i_row)
        try {
            tok_.next_line();
            id = boost::lexical_cast<unsigned int>(*tok_); ++tok_;

            while ( std::min(*bulk_id_iter, *bdr_id_iter) < (int)id) { // skip initialization of some rows in data if ID is missing
                if (*bulk_id_iter < *bdr_id_iter) ++bulk_id_iter;
                else ++bdr_id_iter;
            }

            if (*bulk_id_iter == (int)id) {
                // bulk
                data_cache.read_ascii_data(tok_, header.n_components, (bulk_id_iter - bulk_el_ids.begin()) );
                ++n_bulk_read;  ++bulk_id_iter;
            } else if (*bdr_id_iter == (int)id) {
            	// boundary
                unsigned int bdr_shift = data_cache.get_boundary_begin();
                data_cache.read_ascii_data(tok_, header.n_components, (bdr_id_iter - bdr_el_ids.begin() + bdr_shift) );
                ++n_bdr_read;  ++bdr_id_iter;
            } else {
                if ( (*bulk_id_iter != imax) | (*bdr_id_iter != imax) )
				    WarningOut().fmt("In file '{}', '$ElementData' section for field '{}', time: {}.\nData ID {} is not in order. Skipping rest of data.\n",
                            tok_.f_name(), header.field_name, header.time, id);
                break;
            }

        } catch (boost::bad_lexical_cast &) {
        	THROW(ExcWrongFormat() << EI_Type("$ElementData line") << EI_TokenizerMsg(tok_.position_msg())
        			<< EI_MeshFile(tok_.f_name()) );
        }
    // possibly skip remaining lines after break
    while (i_row < header.n_entities) tok_.next_line(false), ++i_row;

    LogOut().fmt("time: {}; {} bulk and {} boundary entities of field {} read.\n",
    		header.time, n_bulk_read, n_bdr_read, header.field_name);
}



void GmshMeshReader::make_header_table()
{
	header_table_.clear();
	MeshDataHeader header;
	while ( !tok_.eof() ) {
        if ( tok_.skip_to("$ElementData") ) {
            read_data_header(header);
            HeaderTable::iterator it = header_table_.find(header.field_name);

            if (it == header_table_.end()) {  // field doesn't exists, insert new vector to map
            	std::vector<MeshDataHeader> vec;
            	vec.push_back(header);
            	header_table_[header.field_name]=vec;
            } else if ( header.time <= it->second.back().time ) { // time is in wrong order. can't be add
            	WarningOut().fmt("Wrong time order: field '{}', time '{}', file '{}'. Skipping this '$ElementData' section.\n",
            		header.field_name, header.time, tok_.f_name() );
            } else {  // add new time step
            	it->second.push_back(header);
            }
        }
	}

	tok_.set_position( Tokenizer::Position() );
}



BaseMeshReader::MeshDataHeader & GmshMeshReader::find_header(BaseMeshReader::HeaderQuery &header_query)
{
	// check discretization, only type element_data or undefined is supported
	if (header_query.discretization != OutputTime::DiscreteSpace::ELEM_DATA) {
		if (header_query.discretization != OutputTime::DiscreteSpace::UNDEFINED && header_query.discretization != OutputTime::DiscreteSpace::NATIVE_DATA) {
			WarningOut().fmt(
					"Unsupported discretization for field '{}', time: {} and GMSH format.\nType 'ELEM_DATA' of discretization will be used.\n",
					header_query.field_name, header_query.time);
		}
		header_query.discretization = OutputTime::DiscreteSpace::ELEM_DATA;
	}

	HeaderTable::iterator table_it = header_table_.find(header_query.field_name);

	if (table_it == header_table_.end()) {
		// no data found
        THROW( ExcFieldNameNotFound() << EI_FieldName(header_query.field_name) << EI_MeshFile(tok_.f_name()));
	}

	auto comp = [](double t, const MeshDataHeader &a) {
		return t * (1.0 + 2.0*numeric_limits<double>::epsilon()) < a.time;
	};

	std::vector<MeshDataHeader>::iterator headers_it = std::upper_bound(table_it->second.begin(),
			table_it->second.end(),
			header_query.time,
			comp);

	if (headers_it == table_it->second.begin()) {
		THROW( ExcFieldNameNotFound() << EI_FieldName(header_query.field_name)
				                      << EI_MeshFile(tok_.f_name()) << EI_Time(header_query.time));
	}

	--headers_it;
	return *headers_it;
}
