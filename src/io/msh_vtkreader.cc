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
 * @file    msh_vtkreader.cc
 * @brief
 * @author  dalibor
 */


#include <iostream>
#include <vector>
#include "boost/lexical_cast.hpp"

#include "msh_vtkreader.hh"
#include "system/system.hh"
#include "mesh/bih_tree.hh"

#include "config.h"
#include <zlib.h>


/*******************************************************************
 * Helper methods
 */
template<typename T>
T read_binary_value(std::istream &data_stream)
{
	T val;
	data_stream.read(reinterpret_cast<char *>(&val), sizeof(val));
	return val;
}


uint64_t read_header_type(DataType data_header_type, std::istream &data_stream)
{
	if (data_header_type == DataType::uint64)
		return read_binary_value<uint64_t>(data_stream);
	else if (data_header_type == DataType::uint32)
		return (uint64_t)read_binary_value<unsigned int>(data_stream);
	else {
		ASSERT(false).error("Unsupported header_type!\n"); //should not happen
		return 0;
	}
}


/*******************************************************************
 * implementation of VtkMeshReader
 */
const double VtkMeshReader::point_tolerance = 1E-10;


VtkMeshReader::VtkMeshReader(const FilePath &file_name)
: BaseMeshReader(file_name)
{
    data_section_name_ = "DataArray";
    has_compatible_mesh_ = false;
	make_header_table();
}



VtkMeshReader::~VtkMeshReader()
{
	delete data_stream_;
}



void VtkMeshReader::read_base_vtk_attributes(pugi::xml_node vtk_node, unsigned int &n_nodes, unsigned int &n_elements)
{
    try {
    	// header type of appended data
    	header_type_ = this->get_data_type( vtk_node.attribute("header_type").as_string() );
    } catch(ExcWrongType &e) {
    	e << EI_SectionTypeName("base parameter header_type");
    }
	// data format
	if (header_type_ == DataType::undefined) {
		data_format_ = DataFormat::ascii;
	} else {
		if (header_type_!=DataType::uint64 && header_type_!=DataType::uint32) {
			// Allowable values of header type are only 'UInt64' or 'UInt32'
			THROW( ExcWrongType() << EI_ErrMessage("Forbidden") << EI_SectionTypeName("base parameter header_type")
					<< EI_VTKFile(tok_.f_name()));
		}
		std::string compressor = vtk_node.attribute("compressor").as_string();
		if (compressor == "vtkZLibDataCompressor")
			data_format_ = DataFormat::binary_zlib;
		else
			data_format_ = DataFormat::binary_uncompressed;
	}
	// size of node and element vectors
	pugi::xml_node piece_node = vtk_node.child("UnstructuredGrid").child("Piece");
	n_nodes = piece_node.attribute("NumberOfPoints").as_uint();
	n_elements = piece_node.attribute("NumberOfCells").as_uint();
}



Tokenizer::Position VtkMeshReader::get_appended_position() {
	Tokenizer::Position appended_pos;

	{
		// find line by tokenizer
		tok_.set_position( Tokenizer::Position() );
		if (! tok_.skip_to("AppendedData"))
			THROW(ExcMissingTag() << EI_TagType("tag") << EI_TagName("AppendedData") << EI_VTKFile(tok_.f_name()) );
		else {
			appended_pos = tok_.get_position();
		}
	}

	// find exact position of appended data (starts with '_')
	char c;
	data_stream_->seekg(appended_pos.file_position_);
	do {
		data_stream_->get(c);
	} while (c!='_');
	appended_pos.file_position_ = data_stream_->tellg();
	appended_pos.line_counter_++;
	delete data_stream_; // close stream

	// reopen stream in binary mode
	data_stream_ = new std::ifstream( tok_.f_name(), std::ios_base::in | std::ios_base::binary );

	return appended_pos;
}


MeshDataHeader VtkMeshReader::create_header(pugi::xml_node node, unsigned int n_entities, Tokenizer::Position pos)
{
	MeshDataHeader header;
    header.field_name = node.attribute("Name").as_string();
    header.time = 0.0;
    try {
        header.type = this->get_data_type( node.attribute("type").as_string() );
    } catch(ExcWrongType &e) {
        e << EI_SectionTypeName("DataArray " + header.field_name);
    }
    header.n_components = node.attribute("NumberOfComponents").as_uint(1);
    header.n_entities = n_entities;
    std::string format = node.attribute("format").as_string();
    if (format=="appended") {
    	if (data_format_ == DataFormat::ascii)
    		THROW(ExcInvalidFormat() << EI_FieldName(header.field_name) << EI_ExpectedFormat("ascii") << EI_VTKFile(tok_.f_name()) );
        std::streampos file_pos = pos.file_position_;
        file_pos += node.attribute("offset").as_uint();
        header.position = Tokenizer::Position( file_pos, pos.line_counter_, pos.line_position_ );
    } else if (format=="ascii") {
    	if (data_format_ != DataFormat::ascii)
    		THROW(ExcInvalidFormat() << EI_FieldName(header.field_name) << EI_ExpectedFormat("appended") << EI_VTKFile(tok_.f_name()) );

    	tok_.set_position( Tokenizer::Position() );
    	bool is_point = (header.field_name=="");
    	std::string found_str = (is_point) ? "<Points>" : "Name=\"" + header.field_name + "\"";
		if (! tok_.skip_to(found_str))
			THROW(ExcMissingTag() << EI_TagType("DataArray tag") << EI_TagName(header.field_name) << EI_VTKFile(tok_.f_name()) );
		else {
			if (is_point) tok_.skip_to("DataArray");
			header.position = tok_.get_position();
		}
    } else {
    	THROW(ExcUnknownFormat() << EI_FieldName(header.field_name) << EI_VTKFile(tok_.f_name()) );
    }

    return header;
}


void VtkMeshReader::make_header_table()
{
	pugi::xml_document doc;
	doc.load_file( tok_.f_name().c_str() );
	unsigned int n_nodes, n_elements;
	this->read_base_vtk_attributes( doc.child("VTKFile"), n_nodes, n_elements );

	// open ifstream for find position
	data_stream_ = new std::ifstream( tok_.f_name() );

	// data of appended tag
	Tokenizer::Position appended_pos;
	if (header_type_==DataType::undefined) {
		// no AppendedData tag
	} else {
		appended_pos = get_appended_position();
	}

	pugi::xml_node node = doc.child("VTKFile").child("UnstructuredGrid").child("Piece");

	header_table_.clear();

	// create headers of Points and Cells DataArrays
	header_table_["Points"] = create_header( node.child("Points").child("DataArray"), n_nodes, appended_pos );
	header_table_["Points"].field_name = "Points";
	header_table_["connectivity"]
			= create_header( node.child("Cells").find_child_by_attribute("DataArray", "Name", "connectivity"), n_elements, appended_pos );
	header_table_["offsets"]
			= create_header( node.child("Cells").find_child_by_attribute("DataArray", "Name", "offsets"), n_elements, appended_pos );
	header_table_["types"]
			= create_header( node.child("Cells").find_child_by_attribute("DataArray", "Name", "types"), n_elements, appended_pos );

	node = node.child("CellData");
	for (pugi::xml_node subnode = node.child("DataArray"); subnode; subnode = subnode.next_sibling("DataArray")) {
		auto header = create_header( subnode, n_elements, appended_pos );
		header_table_[header.field_name] = header;
	}
}


MeshDataHeader & VtkMeshReader::find_header(double time, std::string field_name)
{
	HeaderTable::iterator table_it = header_table_.find(field_name);

	if (table_it == header_table_.end()) {
		// no data found
        THROW( ExcFieldNameNotFound() << EI_FieldName(field_name) << EI_MeshFile(tok_.f_name()));
	}

	return table_it->second;
}



DataType VtkMeshReader::get_data_type(std::string type_str) {
    // names of types in DataArray section
	static const std::map<std::string, DataType> types = {
			{"Int8",    DataType::int8},
			{"UInt8",   DataType::uint8},
			{"Int16",   DataType::int16},
			{"UInt16",  DataType::uint16},
			{"Int32",   DataType::int32},
			{"UInt32",  DataType::uint32},
			{"Int64",   DataType::int64},
			{"UInt64",  DataType::uint64},
			{"Float32", DataType::float32},
			{"Float64", DataType::float64},
			{"",        DataType::undefined}
	};

	std::map<std::string, DataType>::const_iterator it = types.find(type_str);
	if (it != types.end()) {
	    return it->second;
    } else {
    	THROW( ExcWrongType() << EI_ErrMessage("Unknown") << EI_VTKFile(tok_.f_name()));
        return DataType::uint32;
    }

}



unsigned int VtkMeshReader::type_value_size(DataType data_type)
{
	static const std::vector<unsigned int> sizes = { 1, 1, 2, 2, 4, 4, 8, 8, 4, 8, 0 };

	return sizes[data_type];
}



void VtkMeshReader::read_element_data(ElementDataCacheBase &data_cache, MeshDataHeader actual_header, unsigned int size_of_cache,
		unsigned int n_components, std::vector<int> const & el_ids) {
    switch (data_format_) {
		case DataFormat::ascii: {
			parse_ascii_data( data_cache, size_of_cache, n_components, actual_header.n_entities, actual_header.position );
			break;
		}
		case DataFormat::binary_uncompressed: {
			ASSERT_PTR(data_stream_).error();
			parse_binary_data( data_cache, size_of_cache, n_components, actual_header.n_entities, actual_header.position,
					actual_header.type );
			break;
		}
		case DataFormat::binary_zlib: {
			ASSERT_PTR(data_stream_).error();
			parse_compressed_data( data_cache, size_of_cache, n_components, actual_header.n_entities, actual_header.position,
					actual_header.type);
			break;
		}
		default: {
			ASSERT(false).error(); // should not happen
			break;
		}
	}

    LogOut().fmt("time: {}; {} entities of field {} read.\n",
    		actual_header.time, n_read_, actual_header.field_name);
}


void VtkMeshReader::parse_ascii_data(ElementDataCacheBase &data_cache, unsigned int size_of_cache, unsigned int n_components,
		unsigned int n_entities, Tokenizer::Position pos)
{
    n_read_ = 0;

	tok_.set_position( pos );
    try {
    	tok_.next_line();
    	for (unsigned int i_row = 0; i_row < n_entities; ++i_row) {
    		data_cache.read_ascii_data(tok_, n_components, vtk_to_gmsh_element_map_[i_row]);
            n_read_++;
    	}
	} catch (boost::bad_lexical_cast &) {
		THROW(ExcWrongFormat() << EI_Type("DataArray tag") << EI_TokenizerMsg(tok_.position_msg())
				<< EI_MeshFile(tok_.f_name()) );
	}
}


void VtkMeshReader::parse_binary_data(ElementDataCacheBase &data_cache, unsigned int size_of_cache, unsigned int n_components,
		unsigned int n_entities, Tokenizer::Position pos, DataType value_type)
{
    n_read_ = 0;

    data_stream_->seekg(pos.file_position_);
	uint64_t data_size = read_header_type(header_type_, *data_stream_) / type_value_size(value_type);
	ASSERT_EQ(size_of_cache*n_components*n_entities, data_size).error();

	for (unsigned int i_row = 0; i_row < n_entities; ++i_row) {
		data_cache.read_binary_data(*data_stream_, n_components, vtk_to_gmsh_element_map_[i_row]);
        n_read_++;
	}
}


void VtkMeshReader::parse_compressed_data(ElementDataCacheBase &data_cache, unsigned int size_of_cache, unsigned int n_components,
		unsigned int n_entities, Tokenizer::Position pos, DataType value_type)
{
	data_stream_->seekg(pos.file_position_);
	uint64_t n_blocks = read_header_type(header_type_, *data_stream_);
	uint64_t u_size = read_header_type(header_type_, *data_stream_);
	uint64_t p_size = read_header_type(header_type_, *data_stream_);

	std::vector<uint64_t> block_sizes;
	block_sizes.reserve(n_blocks);
	for (uint64_t i = 0; i < n_blocks; ++i) {
		block_sizes.push_back( read_header_type(header_type_, *data_stream_) );
	}

	stringstream decompressed_data;
	uint64_t decompressed_data_size = 0;
	for (uint64_t i = 0; i < n_blocks; ++i) {
		uint64_t decompressed_block_size = (i==n_blocks-1 && p_size>0) ? p_size : u_size;
		uint64_t compressed_block_size = block_sizes[i];

		std::vector<char> data_block(compressed_block_size);
		data_stream_->read(&data_block[0], compressed_block_size);

		std::vector<char> buffer(decompressed_block_size);

		// set zlib object
		z_stream strm;
		strm.zalloc = 0;
		strm.zfree = 0;
		strm.next_in = (Bytef *)(&data_block[0]);
		strm.avail_in = compressed_block_size;
		strm.next_out = (Bytef *)(&buffer[0]);
		strm.avail_out = decompressed_block_size;

		// decompression of data
		inflateInit(&strm);
		inflate(&strm, Z_NO_FLUSH);
		inflateEnd(&strm);

		// store decompressed data to stream
		decompressed_data << std::string(buffer.begin(), buffer.end());
		decompressed_data_size += decompressed_block_size;
	}

    n_read_ = 0;

	uint64_t data_size = decompressed_data_size / type_value_size(value_type);
	ASSERT_EQ(size_of_cache*n_components*n_entities, data_size).error();

	for (unsigned int i_row = 0; i_row < n_entities; ++i_row) {
		data_cache.read_binary_data(decompressed_data, n_components, vtk_to_gmsh_element_map_[i_row]);
        n_read_++;
	}
}


void VtkMeshReader::check_compatible_mesh(Mesh &mesh)
{
    std::vector<int> el_ids;
    std::vector<unsigned int> node_ids; // allow mapping ids of nodes from VTK mesh to GMSH
    std::vector<unsigned int> offsets_vec; // value of offset section in VTK file

    {
        // read points data section, find corresponding nodes in GMSH trough BIH tree
        // points in data section and nodes in GMSH must be in ratio 1:1
        // store orders (mapping between VTK and GMSH file) into node_ids vector
        std::vector<unsigned int> searched_elements; // for BIH tree
        unsigned int i_node, i_elm_node;
        const BIHTree &bih_tree=mesh.get_bih_tree();

        NodeDataTable nodes_data = this->read_nodes_data();
        node_ids.resize(nodes_data.size());

        for (node_data : nodes_data) {
            int found_i_node = -1;
            bih_tree.find_point(node_data.second, searched_elements);

            for (std::vector<unsigned int>::iterator it = searched_elements.begin(); it!=searched_elements.end(); it++) {
                ElementFullIter ele = mesh.element( *it );
                FOR_ELEMENT_NODES(ele, i_node)
                {
                    if ( compare_points(ele->node[i_node]->point(), node_data.second) ) {
                    	i_elm_node = mesh.node_vector.index(ele->node[i_node]);
                        if (found_i_node == -1) found_i_node = i_elm_node;
                        else if (found_i_node != i_elm_node) {
                        	THROW( ExcIncompatibleMesh() << EI_ErrMessage("duplicate nodes found in GMSH file")
                        			<< EI_VTKFile(tok_.f_name()));
                        }
                    }
                }
            }
            if (found_i_node == -1) {
            	THROW( ExcIncompatibleMesh() << EI_ErrMessage("no node found in GMSH file") << EI_VTKFile(tok_.f_name()));
            }
            node_ids[ node_data.first ] = (unsigned int)found_i_node;
            searched_elements.clear();
        }

    }

    {
        // read offset data section into offsets_vec vector, it's used for reading connectivity
        MeshDataHeader offset_header = this->find_header(0.0, "offsets");
        for (unsigned int i=el_ids.size(); i<offset_header.n_entities; ++i) {
        	el_ids.push_back(i);
        	vtk_to_gmsh_element_map_.push_back(i);
        }

        ElementDataCache<unsigned int> offset_cache(offset_header, 1, offset_header.n_components*offset_header.n_entities);
        this->read_element_data(offset_cache, offset_header, 1, offset_header.n_components, el_ids);

        offsets_vec = *(offset_cache.get_component_data(0) );
    }

    {
        // read connectivity data section, find corresponding elements in GMSH
        // cells in data section and elements in GMSH must be in ratio 1:1
        // store orders (mapping between VTK and GMSH file) into vtk_to_gmsh_element_map_ vector
        MeshDataHeader con_header = this->find_header(0.0, "connectivity");
        con_header.n_entities = offsets_vec[offsets_vec.size()-1];
        for (unsigned int i=el_ids.size(); i<con_header.n_entities; ++i) {
        	el_ids.push_back(i);
        	vtk_to_gmsh_element_map_.push_back(i);
        }

        ElementDataCache<unsigned int> con_cache(con_header, 1, con_header.n_components*con_header.n_entities);
        this->read_element_data(con_cache, con_header, 1, con_header.n_components, el_ids);

        std::vector<unsigned int> &connectivity_vec = *(con_cache.get_component_data(0) );
        vector<unsigned int> node_list;
        vector<unsigned int> candidate_list; // returned by intersect_element_lists
        vector<unsigned int> result_list; // list of elements with same dimension as vtk element
        vtk_to_gmsh_element_map_.clear();
        vtk_to_gmsh_element_map_.resize(offsets_vec.size());
        // iterate trough connectivity data, to each VTK cell must exist only one GMSH element
        // fill vtk_to_gmsh_element_map_ vector
        unsigned int i_con = 0, last_offset=0, dim;
        for (unsigned int i=0; i<offsets_vec.size(); ++i) { // iterate trough offset - one value for every element
        	dim = offsets_vec[i] - last_offset - 1; // dimension of vtk element
            for ( ; i_con<offsets_vec[i]; ++i_con ) { // iterate trough all nodes of any element
                node_list.push_back( node_ids[connectivity_vec[i_con]] );
            }
            mesh.intersect_element_lists(node_list, candidate_list);
            for (auto i_elm : candidate_list) {
            	if ( mesh.element( i_elm )->dim() == dim ) result_list.push_back(i_elm);
            }
            if (result_list.size() != 1) {
            	THROW( ExcIncompatibleMesh() << EI_ErrMessage("intersect_element_lists must produce one element")
            			<< EI_VTKFile(tok_.f_name()));
            }
            vtk_to_gmsh_element_map_[i] = result_list[0];
            node_list.clear();
            result_list.clear();
            last_offset = offsets_vec[i];
        }
    }

    has_compatible_mesh_ = true;
}


NodeDataTable VtkMeshReader::read_nodes_data() {
    using namespace boost;
    NodeDataTable node_data_table;
    std::vector<int> el_ids;

    MeshDataHeader point_header = this->find_header(0.0, "Points");
    ASSERT_EQ(3, point_header.n_components).error();
    // fill vectors necessary for correct reading of data, we read all data same order as data is stored
    for (unsigned int i=0; i<point_header.n_entities; ++i) {
    	el_ids.push_back(i);
    	vtk_to_gmsh_element_map_.push_back(i);
    }

    // create temporary data cache
    ElementDataCache<double> node_cache(point_header, 1, point_header.n_components*point_header.n_entities);

    // check compatible nodes, to each VTK point must exist only one GMSH node
    this->read_element_data(node_cache, point_header, 1, point_header.n_components, el_ids);
    std::vector<double> &node_vec = *(node_cache.get_component_data(0) );
    ASSERT_EQ(node_vec.size(), point_header.n_components*point_header.n_entities).error();
    for (unsigned int i=0; i<point_header.n_entities; ++i) { // node id is equivalent with 'i'
    	arma::vec3 node_data;                                // node coordinates
    	node_data(0) = node_vec[3*i];
    	node_data(1) = node_vec[3*i+1];
    	node_data(2) = node_vec[3*i+2];
        node_data_table.push_back( std::make_pair(i, node_data) );
    }
    vtk_to_gmsh_element_map_.clear();

    return node_data_table;
}


bool VtkMeshReader::compare_points(arma::vec3 &p1, arma::vec3 &p2) {
	return fabs(p1[0]-p2[0]) < point_tolerance
		&& fabs(p1[1]-p2[1]) < point_tolerance
		&& fabs(p1[2]-p2[2]) < point_tolerance;
}



// explicit instantiation of template methods
template unsigned int read_binary_value<unsigned int>(std::istream &data_stream);
template uint64_t read_binary_value<uint64_t>(std::istream &data_stream);
