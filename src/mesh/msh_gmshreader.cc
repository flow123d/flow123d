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

#include <istream>
#include <string>
#include <limits>

#include "msh_gmshreader.h"

#include "system/system.hh"
#include "system/sys_profiler.hh"
#include "system/tokenizer.hh"
#include "boost/lexical_cast.hpp"

#include "mesh/mesh.h"
#include "mesh/nodes.hh"


using namespace std;


GmshMeshReader::GmshMeshReader(const FilePath &file_name)
: tok_(file_name)
{
	current_cache_ = new ElementDataCacheBase();
    tok_.set_comment_pattern( "#");
    make_header_table();
}



GmshMeshReader::GmshMeshReader(std::istream &in)
: tok_(in)
{
	current_cache_ = new ElementDataCacheBase();
    tok_.set_comment_pattern( "#");
    make_header_table();
}



GmshMeshReader::~GmshMeshReader()   // Tokenizer close the file automatically
{}



void GmshMeshReader::read_mesh(Mesh* mesh, const RegionDB::MapElementIDToRegionID *el_to_reg_map) {
    START_TIMER("GMSHReader - read mesh");
    
    ASSERT( mesh , "Argument mesh is NULL.\n");
    read_physical_names(tok_, mesh);
    read_nodes(tok_, mesh);
    read_elements(tok_, mesh, el_to_reg_map);
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
        for (unsigned int i = 0; i < n_nodes; ++i) {
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



void GmshMeshReader::read_elements(Tokenizer &tok, Mesh * mesh, const RegionDB::MapElementIDToRegionID *el_to_reg_map) {
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
            // 15 Point (1 node)
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
                case 15:
                    dim = 0;
                    break;
                default:
                    dim = 0;
                    xprintf(UsrErr, "Element %d is of the unsupported type %d\n", id, type);
                    break;
            }

            //get number of tags (at least 2)
            unsigned int n_tags = lexical_cast<unsigned int>(*tok);
            INPUT_CHECK(n_tags >= 2, "At least two element tags have to be defined for element with id=%d, %s.\n",
                    id, tok.position_msg().c_str());
            ++tok;

            //get tags 1 and 2
            unsigned int region_id = lexical_cast<unsigned int>(*tok); ++tok;
            lexical_cast<unsigned int>(*tok); ++tok; // GMSH region number, we do not store this
            //get remaining tags
            unsigned int partition_id=0;
            if (n_tags > 2)  { partition_id = lexical_cast<unsigned int>(*tok); ++tok; } // save partition number from the new GMSH format
            for (unsigned int ti = 3; ti < n_tags; ti++) ++tok;         //skip remaining tags

            // allocate element arrays TODO: should be in mesh class
            Element *ele=nullptr;
            // possibly modify region id
            if (el_to_reg_map) {
                RegionDB::MapElementIDToRegionID::const_iterator it = el_to_reg_map->find(id);
                if (it != el_to_reg_map->end()) region_id = it->second;
            }
            RegionIdx region_idx = mesh->region_db_.add_region( region_id, dim );

            if (region_idx.is_boundary()) {
                ele = mesh->bc_elements.add_item(id);
            } else {
                if(dim == 0 )
                    xprintf(Warn, "Bulk elements of zero size(dim=0) are not supported. Mesh file: %s, Element ID: %d.\n", tok.f_name().c_str() ,id);
                else
                    ele = mesh->element.add_item(id);
            }
            ele->init(dim, mesh, region_idx);
            ele->pid=partition_id;

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

    mesh->n_all_input_elements_=mesh->element.size() + mesh->bc_elements.size();
    xprintf(Msg, " %d bulk elements, %d boundary elements. \n", mesh->element.size(), mesh->bc_elements.size());
}



void GmshMeshReader::read_physical_names(Tokenizer &tok, Mesh * mesh) {
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

            mesh->region_db_.add_region(id, name, dim);
        }

    } catch (bad_lexical_cast &) {
        xprintf(UsrErr, "Wrong format of number, %s.\n", tok.position_msg().c_str());
    }
}


// Is assumed to be called just after tok.skip_to("..")
// reads the header from the tokenizer @p tok and return it as the second parameter
void GmshMeshReader::read_data_header(Tokenizer &tok, GMSH_DataHeader &head) {
    using namespace boost;
    try {
        // string tags
        tok.next_line(false);
        unsigned int n_str = lexical_cast<unsigned int>(*tok); ++tok;
        head.field_name="";
        head.interpolation_scheme = "";
        if (n_str > 0) {
            tok.next_line(); n_str--;
            head.field_name= *tok; ++tok; //  unquoted by tokenizer if needed
        }
        if (n_str > 0) {
            tok.next_line(); n_str--;
            head.interpolation_scheme = *tok; ++tok;
        }
        for(;n_str>0;n_str--) tok.next_line(false); // skip possible remaining tags

        //real tags
        tok.next_line();
        unsigned int n_real = lexical_cast<unsigned int>(*tok); ++tok;
        head.time=0.0;
        if (n_real>0) {
            tok.next_line(); n_real--;
            head.time=lexical_cast<double>(*tok); ++tok;
        }
        for(;n_real>0;n_real--) tok.next_line(false);

        // int tags
        tok.next_line();
        unsigned int n_int = lexical_cast<unsigned int>(*tok); ++tok;
        head.time_index=0;
        head.n_components=1;
        head.n_entities=0;
        head.partition_index=0;
        if (n_int>0) {
            tok.next_line(); n_int--;
            head.time_index=lexical_cast<unsigned int>(*tok); ++tok;
        }
        if (n_int>0) {
            tok.next_line(); n_int--;
            head.n_components=lexical_cast<unsigned int>(*tok); ++tok;
        }
        if (n_int>0) {
            tok.next_line(); n_int--;
            head.n_entities=lexical_cast<unsigned int>(*tok); ++tok;
        }
        for(;n_int>0;n_int--) tok.next_line(false);
        head.position = tok_.get_position();
    } catch (bad_lexical_cast &) {
                xprintf(UsrErr, "Wrong format of the $ElementData header, %s.\n", tok.position_msg().c_str());
    }
}



void GmshMeshReader::read_element_data( GMSH_DataHeader &search_header,
        double *data, std::vector<int> const & el_ids)
{

    using namespace boost;

    unsigned int id, idx, i_row;
    unsigned int n_read = 0;
    vector<int>::const_iterator id_iter = el_ids.begin();
    double * data_ptr;
    GMSH_DataHeader actual_header = find_header(search_header.time, search_header.field_name);

    // check that the header is valid, try to correct
    if (actual_header.n_components != search_header.n_components) {
        xprintf(Warn, "In file '%s', '$ElementData' section for field '%s', time: %f.\nWrong number of components: %d, using %d instead.\n",
                tok_.f_name().c_str(), search_header.field_name.c_str(), actual_header.time, actual_header.n_components, search_header.n_components);
        actual_header.n_components=search_header.n_components;
    }
    if (actual_header.n_entities != search_header.n_entities) {
        xprintf(Warn, "In file '%s', '$ElementData' section for field '%s', time: %f.\nWrong number of entities: %d, using %d instead.\n",
                tok_.f_name().c_str(), search_header.field_name.c_str(), actual_header.time, actual_header.n_entities, search_header.n_entities);
        // actual_header.n_entities=search_header.n_entities;
    }

    // read @p data buffer as we have correct header with already passed time
    // we assume that @p data buffer is big enough
    tok_.set_position(actual_header.position);
    for (i_row = 0; i_row < actual_header.n_entities; ++i_row)
        try {
            tok_.next_line();
            id = lexical_cast<unsigned int>(*tok_); ++tok_;
            while (id_iter != el_ids.end() && *id_iter < (int)id) {
                ++id_iter; // skip initialization of some rows in data if ID is missing
            }
            if (id_iter == el_ids.end()) {
                xprintf(Warn,"In file '%s', '$ElementData' section for field '%s', time: %f.\nData ID %d not found or is not in order. Skipping rest of data.\n",
                        tok_.f_name().c_str(), search_header.field_name.c_str(), actual_header.time, id);
                break;
            }
            // save data from the line if ID was found
            if (*id_iter == (int)id) {
                idx = id_iter - el_ids.begin();
                data_ptr = data + idx * search_header.n_components;
                for (unsigned int i_col =0; i_col < search_header.n_components; ++i_col, ++data_ptr) {
                    *(data_ptr) = lexical_cast<double>(*tok_); ++tok_;
                }
                n_read++;
            }
            // skip the line if ID on the line  < actual ID in the map el_ids
        } catch (bad_lexical_cast &) {
            xprintf(UsrErr, "Wrong format of $ElementData line, %s.\n", tok_.position_msg().c_str());
        }
    // possibly skip remaining lines after break
    while (i_row < actual_header.n_entities) tok_.next_line(false), ++i_row;

    xprintf(MsgLog, "time: %f; %d entities of field %s read.\n",
    		actual_header.time, n_read, actual_header.field_name.c_str());

    search_header.actual = true; // use input header to indicate modification of @p data buffer
}



template<typename T>
typename ElementDataCache<T>::ComponentDataPtr GmshMeshReader::get_element_data( GMSH_DataHeader &search_header,
		std::vector<int> const & el_ids, unsigned int component_idx)
{
    using namespace boost;

    if ( !current_cache_->is_actual(search_header.time, search_header.field_name) ) {

	    unsigned int id, idx, i_row;
	    unsigned int n_read = 0;
    	unsigned int size_of_cache; // count of vectors stored in cache
	    vector<int>::const_iterator id_iter = el_ids.begin();
	    GMSH_DataHeader actual_header = find_header(search_header.time, search_header.field_name);

	    // check that the header is valid, try to correct
	    if (actual_header.n_entities != search_header.n_entities) {
	        xprintf(Warn, "In file '%s', '$ElementData' section for field '%s', time: %f.\nWrong number of entities: %d, using %d instead.\n",
	                tok_.f_name().c_str(), search_header.field_name.c_str(), actual_header.time, actual_header.n_entities, search_header.n_entities);
	        // actual_header.n_entities=search_header.n_entities;
	    }

	    if (search_header.n_components == 1) {
	    	// read for MultiField to 'n_comp' vectors
	    	// or for Field if ElementData contains only one value
	    	size_of_cache = actual_header.n_components;
	    }
	    else {
	    	// read for Field if more values is stored to one vector
	    	size_of_cache = 1;
	    	if (actual_header.n_components != search_header.n_components) {
		        xprintf(Warn, "In file '%s', '$ElementData' section for field '%s', time: %f.\nWrong number of components: %d, using %d instead.\n",
		                tok_.f_name().c_str(), search_header.field_name.c_str(), actual_header.time, actual_header.n_components, search_header.n_components);
		        actual_header.n_components=search_header.n_components;
	    	}
	    }

	    // create vector of shared_ptr for cache
	    typename ElementDataCache<T>::CacheData data_cache(size_of_cache);
	    for (unsigned int i=0; i<size_of_cache; ++i) {
			typename ElementDataCache<T>::ComponentDataPtr row_vec = std::make_shared<std::vector<T>>();
			row_vec->resize(search_header.n_components*search_header.n_entities);
			data_cache[i] = row_vec;
	    }

	    // read @p data buffer as we have correct header with already passed time
	    // we assume that @p data buffer is big enough
	    tok_.set_position(actual_header.position);

	    // read data
	    for (i_row = 0; i_row < actual_header.n_entities; ++i_row)
	        try {
	            tok_.next_line();
	            id = lexical_cast<unsigned int>(*tok_); ++tok_;
	            //skip_element = false;
	            while (id_iter != el_ids.end() && *id_iter < (int)id) {
	                ++id_iter; // skip initialization of some rows in data if ID is missing
	            }
	            if (id_iter == el_ids.end()) {
	                xprintf(Warn,"In file '%s', '$ElementData' section for field '%s', time: %f.\nData ID %d not found or is not in order. Skipping rest of data.\n",
	                        tok_.f_name().c_str(), search_header.field_name.c_str(), actual_header.time, id);
	                break;
	            }
	            // save data from the line if ID was found
	            if (*id_iter == (int)id) {
	            	for (unsigned int i_vec=0; i_vec<size_of_cache; ++i_vec) {
	            		idx = (id_iter - el_ids.begin()) * search_header.n_components;
	            		std::vector<T> &vec = *( data_cache[i_vec].get() );
	            		for (unsigned int i_col=0; i_col < search_header.n_components; ++i_col, ++idx) {
	            			vec[idx] = lexical_cast<T>(*tok_);
	            			++tok_;
	            		}
	            	}
	                n_read++;
	            }
	            // skip the line if ID on the line  < actual ID in the map el_ids
	        } catch (bad_lexical_cast &) {
	            xprintf(UsrErr, "Wrong format of $ElementData line, %s.\n", tok_.position_msg().c_str());
	        }
	    // possibly skip remaining lines after break
	    while (i_row < actual_header.n_entities) tok_.next_line(false), ++i_row;

	    xprintf(Msg, "time: %f; %d entities of field %s read.\n",
	    		actual_header.time, n_read, actual_header.field_name.c_str());

	    search_header.actual = true; // use input header to indicate modification of @p data buffer

	    // set new cache
	    delete current_cache_;
	    current_cache_ = new ElementDataCache<T>(actual_header.time, actual_header.field_name, data_cache);
	}

    if (component_idx == std::numeric_limits<unsigned int>::max()) component_idx = 0;
	return static_cast< ElementDataCache<T> *>(current_cache_)->get_component_data(component_idx);
}



void GmshMeshReader::make_header_table()
{
	header_table_.clear();
	GMSH_DataHeader header;
	while ( !tok_.eof() ) {
        if ( tok_.skip_to("$ElementData") ) {
            read_data_header(tok_, header);
            HeaderTable::iterator it = header_table_.find(header.field_name);

            if (it == header_table_.end()) {  // field doesn't exists, insert new vector to map
            	std::vector<GMSH_DataHeader> vec;
            	vec.push_back(header);
            	//header_table_.insert( std::pair<std::string, std::vector<GMSH_DataHeader> >(header.field_name, vec) );
            	header_table_[header.field_name]=vec;
            } else if ( header.time <= it->second.back().time ) {  // time is in wrong order. can't be add
            	xprintf(Warn,
            		"Non-ascending time series detected in file '%s', '$ElementData' section, quantity: '%s', time: '%d'. Skipping this time.\n",
            		tok_.f_name().c_str(), header.field_name.c_str(), header.time);
            } else {  // add new time step
            	it->second.push_back(header);
            }
        }
	}

	tok_.set_position( Tokenizer::Position() );
}



GMSH_DataHeader &  GmshMeshReader::find_header(double time, std::string field_name)
{
	HeaderTable::iterator table_it = header_table_.find(field_name);

	if (table_it == header_table_.end()) {
		// no data found
        THROW( ExcFieldNameNotFound() << EI_FieldName(field_name) << EI_GMSHFile(tok_.f_name()));
	}

	auto comp = [](double t, const GMSH_DataHeader &a) {
		return t * (1.0 + 2.0*numeric_limits<double>::epsilon()) < a.time;
	};

	std::vector<GMSH_DataHeader>::iterator headers_it = std::upper_bound(table_it->second.begin(),
			table_it->second.end(),
			time,
			comp);

	if (headers_it == table_it->second.begin()) {
		THROW( ExcFieldNameNotFound() << EI_FieldName(field_name)
				                      << EI_GMSHFile(tok_.f_name()) << EI_Time(time));
	}

	--headers_it;
	return *headers_it;
}


// explicit instantiation of template methods
#define READER_GET_ELEMENT_DATA(TYPE) \
template typename ElementDataCache<TYPE>::ComponentDataPtr GmshMeshReader::get_element_data<TYPE>(GMSH_DataHeader &search_header, \
	std::vector<int> const & el_ids, unsigned int component_idx)

READER_GET_ELEMENT_DATA(int);
READER_GET_ELEMENT_DATA(unsigned int);
READER_GET_ELEMENT_DATA(double);
