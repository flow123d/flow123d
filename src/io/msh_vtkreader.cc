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
#include <pugixml.hpp>
#include "boost/lexical_cast.hpp"

#include "msh_vtkreader.hh"
#include "system/system.hh"
#include "system/index_types.hh"
#include "mesh/bih_tree.hh"
#include "mesh/mesh.h"
#include "mesh/accessors.hh"

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
: BaseMeshReader(file_name),
  time_step_(0.0)
{
    data_section_name_ = "DataArray";
    has_compatible_mesh_ = false;
	make_header_table();
}



VtkMeshReader::VtkMeshReader(const FilePath &file_name, std::shared_ptr<ElementDataFieldMap> element_data_values, double time_step)
: BaseMeshReader(file_name, element_data_values),
  time_step_(time_step)
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


BaseMeshReader::MeshDataHeader VtkMeshReader::create_header(pugi::xml_node node, unsigned int n_entities, Tokenizer::Position pos,
		OutputTime::DiscreteSpace disc)
{
	MeshDataHeader header;
    header.field_name = node.attribute("Name").as_string();
    header.time = this->time_step_;
    try {
        header.type = this->get_data_type( node.attribute("type").as_string() );
    } catch(ExcWrongType &e) {
        e << EI_SectionTypeName("DataArray " + header.field_name);
    }
    header.discretization = disc;
    if (disc == OutputTime::DiscreteSpace::NATIVE_DATA) {
        header.n_components = node.attribute("n_dofs_per_element").as_uint();
    	header.dof_handler_hash = node.attribute("dof_handler_hash").as_uint();
    } else {
        header.n_components = node.attribute("NumberOfComponents").as_uint(1);
    }
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
	auto points_header = create_header( node.child("Points").child("DataArray"), n_nodes, appended_pos,
			OutputTime::DiscreteSpace::MESH_DEFINITION );
	points_header.field_name = "Points";
	header_table_.insert( std::pair<std::string, MeshDataHeader>("Points", points_header) );
	auto con_header = create_header( node.child("Cells").find_child_by_attribute("DataArray", "Name", "connectivity"),
			n_elements, appended_pos, OutputTime::DiscreteSpace::MESH_DEFINITION );
	header_table_.insert( std::pair<std::string, MeshDataHeader>("connectivity", con_header) );
	auto offsets_header = create_header( node.child("Cells").find_child_by_attribute("DataArray", "Name", "offsets"),
			n_elements, appended_pos, OutputTime::DiscreteSpace::MESH_DEFINITION );
	header_table_.insert( std::pair<std::string, MeshDataHeader>("offsets", offsets_header) );
	auto types_header = create_header( node.child("Cells").find_child_by_attribute("DataArray", "Name", "types"), n_elements,
			appended_pos, OutputTime::DiscreteSpace::MESH_DEFINITION );
	header_table_.insert( std::pair<std::string, MeshDataHeader>("types", types_header) );

	pugi::xml_node point_node = node.child("PointData");
	for (pugi::xml_node subnode = point_node.child("DataArray"); subnode; subnode = subnode.next_sibling("DataArray")) {
		auto header = create_header( subnode, n_nodes, appended_pos, OutputTime::DiscreteSpace::NODE_DATA );
		header_table_.insert( std::pair<std::string, MeshDataHeader>(header.field_name, header) );
	}

	pugi::xml_node cell_node = node.child("CellData");
	for (pugi::xml_node subnode = cell_node.child("DataArray"); subnode; subnode = subnode.next_sibling("DataArray")) {
		auto header = create_header( subnode, n_elements, appended_pos, OutputTime::DiscreteSpace::ELEM_DATA );
		header_table_.insert( std::pair<std::string, MeshDataHeader>(header.field_name, header) );
	}

	pugi::xml_node native_node = node.child("Flow123dData");
	for (pugi::xml_node subnode = native_node.child("DataArray"); subnode; subnode = subnode.next_sibling("DataArray")) {
		auto header = create_header( subnode, n_elements, appended_pos, OutputTime::DiscreteSpace::NATIVE_DATA );
		header_table_.insert( std::pair<std::string, MeshDataHeader>(header.field_name, header) );
	}
}


BaseMeshReader::MeshDataHeader & VtkMeshReader::find_header(BaseMeshReader::HeaderQuery &header_query)
{
	unsigned int count = header_table_.count(header_query.field_name);

	if (count == 0) {
		// no data found
        THROW( ExcFieldNameNotFound() << EI_FieldName(header_query.field_name) << EI_MeshFile(tok_.f_name()));
	} else if (count == 1) {
		HeaderTable::iterator table_it = header_table_.find(header_query.field_name);

		// check discretization
		if (header_query.discretization != table_it->second.discretization) {
			if (header_query.discretization != OutputTime::DiscreteSpace::UNDEFINED) {
				WarningOut().fmt(
						"Invalid value of 'input_discretization' for field '{}', time: {}.\nCorrect discretization type will be used.\n",
						header_query.field_name, header_query.time);
			}
			header_query.discretization = table_it->second.discretization;
		}

		if (header_query.discretization == OutputTime::DiscreteSpace::NATIVE_DATA)
			if (header_query.dof_handler_hash != table_it->second.dof_handler_hash) {
				THROW(ExcInvalidDofHandler() << EI_FieldName(header_query.field_name) << EI_VTKFile(tok_.f_name()) );
			}
		actual_header_ = table_it->second;
	} else {
		/*HeaderTable::iterator table_it;
		for (table_it=header_table_.equal_range(header_query.field_name).first; table_it!=header_table_.equal_range(header_query.field_name).second; ++table_it) {
			if (header_query.discretization != table_it->second.discretization) {
				header_query.dof_handler_hash = table_it->second.dof_handler_hash;
				actual_header_ = table_it->second;
			}
		}*/
		THROW( ExcMissingFieldDiscretization() << EI_FieldName(header_query.field_name) << EI_Time(header_query.time) << EI_MeshFile(tok_.f_name()));
	}

	return actual_header_;
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



void VtkMeshReader::read_element_data(ElementDataCacheBase &data_cache, MeshDataHeader actual_header, unsigned int n_components,
		bool boundary_domain) {

	ASSERT(!boundary_domain).error("Reading VTK data of boundary elements is not supported yet!\n");

    switch (data_format_) {
		case DataFormat::ascii: {
			parse_ascii_data( data_cache, n_components, actual_header.n_entities, actual_header.position, boundary_domain );
			break;
		}
		case DataFormat::binary_uncompressed: {
			ASSERT_PTR(data_stream_).error();
			parse_binary_data( data_cache, n_components, actual_header.n_entities, actual_header.position, boundary_domain);
			break;
		}
		case DataFormat::binary_zlib: {
			ASSERT_PTR(data_stream_).error();
			parse_compressed_data( data_cache, n_components, actual_header.n_entities, actual_header.position, boundary_domain);
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


void VtkMeshReader::parse_ascii_data(ElementDataCacheBase &data_cache, unsigned int n_components, unsigned int n_entities,
		Tokenizer::Position pos, bool boundary_domain)
{
    n_read_ = 0;

	tok_.set_position( pos );
    try {
    	tok_.next_line();
    	for (unsigned int i_row = 0; i_row < n_entities; ++i_row) {
    		data_cache.read_ascii_data(tok_, n_components, get_element_vector(boundary_domain)[i_row]);
            n_read_++;
    	}
	} catch (boost::bad_lexical_cast &) {
		THROW(ExcWrongFormat() << EI_Type("DataArray tag") << EI_TokenizerMsg(tok_.position_msg())
				<< EI_MeshFile(tok_.f_name()) );
	}
}


void VtkMeshReader::parse_binary_data(ElementDataCacheBase &data_cache, unsigned int n_components, unsigned int n_entities,
		Tokenizer::Position pos, bool boundary_domain)
{
    n_read_ = 0;

    data_stream_->seekg(pos.file_position_);
	read_header_type(header_type_, *data_stream_);

	for (unsigned int i_row = 0; i_row < n_entities; ++i_row) {
		data_cache.read_binary_data(*data_stream_, n_components, get_element_vector(boundary_domain)[i_row]);
        n_read_++;
	}
}


void VtkMeshReader::parse_compressed_data(ElementDataCacheBase &data_cache, unsigned int n_components, unsigned int n_entities,
		Tokenizer::Position pos, bool boundary_domain)
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

	for (unsigned int i_row = 0; i_row < n_entities; ++i_row) {
		data_cache.read_binary_data(decompressed_data, n_components, get_element_vector(boundary_domain)[i_row]);
        n_read_++;
	}
}


void VtkMeshReader::check_compatible_mesh(Mesh &mesh)
{
	this->create_node_element_caches();
    std::vector<unsigned int> node_ids; // allow mapping ids of nodes from VTK mesh to GMSH
    std::vector<unsigned int> offsets_vec; // value of offset section in VTK file

	bulk_elements_id_.clear();

    {
        // read points data section, find corresponding nodes in GMSH trough BIH tree
        // points in data section and nodes in GMSH must be in ratio 1:1
        // store orders (mapping between VTK and GMSH file) into node_ids vector
        std::vector<unsigned int> searched_elements; // for BIH tree
        unsigned int i_node, i_elm_node;
        const BIHTree &bih_tree=mesh.get_bih_tree();

    	ElementDataFieldMap::iterator field_it=element_data_values_->find("Points");
    	ASSERT(field_it != element_data_values_->end()).error("Missing cache of Points section. Did you call create_node_element_caches()?\n");

    	// create nodes of mesh
    	std::vector<double> &node_vec = *( dynamic_cast<ElementDataCache<double> &>(*(field_it->second)).get_component_data(0).get() );
    	unsigned int n_nodes = node_vec.size() / 3;
    	node_ids.resize(n_nodes);
        for (unsigned int i=0; i<n_nodes; ++i) {
            arma::vec3 point = { node_vec[3*i], node_vec[3*i+1], node_vec[3*i+2] };
            uint found_i_node = Mesh::undef_idx;
            bih_tree.find_point(point, searched_elements);

            for (std::vector<unsigned int>::iterator it = searched_elements.begin(); it!=searched_elements.end(); it++) {
                ElementAccessor<3> ele = mesh.element_accessor( *it );
                for (i_node=0; i_node<ele->n_nodes(); i_node++)
                {
                    if ( compare_points(*ele.node(i_node), point) ) {
                    	i_elm_node = ele.node(i_node).idx();
                        if (found_i_node == Mesh::undef_idx) found_i_node = i_elm_node;
                        else if (found_i_node != i_elm_node) {
                        	THROW( ExcIncompatibleMesh() << EI_ErrMessage("duplicate nodes found in GMSH file")
                        			<< EI_VTKFile(tok_.f_name()));
                        }
                    }
                }
            }
            if (found_i_node == Mesh::undef_idx) {
            	THROW( ExcIncompatibleMesh() << EI_ErrMessage("no node found in GMSH file") << EI_VTKFile(tok_.f_name()));
            }
            node_ids[i] = found_i_node;
            searched_elements.clear();
        }

    }

    {
    	ElementDataFieldMap::iterator field_it=element_data_values_->find("offsets");
    	ASSERT(field_it != element_data_values_->end()).error("Missing cache of offsets section. Did you call create_node_element_caches()?\n");

        offsets_vec = *( dynamic_cast<ElementDataCache<unsigned int> &>(*(field_it->second)).get_component_data(0).get() );
    }

    {
        // read connectivity data section, find corresponding elements in GMSH
        // cells in data section and elements in GMSH must be in ratio 1:1
        // store orders (mapping between VTK and GMSH file) into bulk_elements_id_ vector
    	ElementDataFieldMap::iterator field_it=element_data_values_->find("connectivity");
    	ASSERT(field_it != element_data_values_->end()).error("Missing cache of connectivity section. Did you call create_node_element_caches()?\n");

        std::vector<unsigned int> &connectivity_vec = *( dynamic_cast<ElementDataCache<unsigned int> &>(*(field_it->second)).get_component_data(0).get() );
        vector<unsigned int> node_list;
        vector<unsigned int> candidate_list; // returned by intersect_element_lists
        vector<unsigned int> result_list; // list of elements with same dimension as vtk element
        bulk_elements_id_.clear();
        bulk_elements_id_.resize(offsets_vec.size());
        // iterate trough connectivity data, to each VTK cell must exist only one GMSH element
        // fill bulk_elements_id_ vector
        unsigned int i_con = 0, last_offset=0, dim;
        for (unsigned int i=0; i<offsets_vec.size(); ++i) { // iterate trough offset - one value for every element
        	dim = offsets_vec[i] - last_offset - 1; // dimension of vtk element
            for ( ; i_con<offsets_vec[i]; ++i_con ) { // iterate trough all nodes of any element
                node_list.push_back( node_ids[connectivity_vec[i_con]] );
            }
            mesh.intersect_element_lists(node_list, candidate_list);
            for (auto i_elm : candidate_list) {
            	if ( mesh.element_accessor(i_elm)->dim() == dim ) result_list.push_back(i_elm);
            }
            if (result_list.size() != 1) {
            	THROW( ExcIncompatibleMesh() << EI_ErrMessage("intersect_element_lists must produce one element")
            			<< EI_VTKFile(tok_.f_name()));
            }
            bulk_elements_id_[i] = (LongIdx)result_list[0];
            node_list.clear();
            result_list.clear();
            last_offset = offsets_vec[i];
        }
    }

    has_compatible_mesh_ = true;
}


bool VtkMeshReader::compare_points(const arma::vec3 &p1, const arma::vec3 &p2) {
	return fabs(p1[0]-p2[0]) < point_tolerance
		&& fabs(p1[1]-p2[1]) < point_tolerance
		&& fabs(p1[2]-p2[2]) < point_tolerance;
}


void VtkMeshReader::read_physical_names(Mesh*) {
	// will be implemented later
	// ASSERT(0).error("Not implemented!");
}


void VtkMeshReader::read_nodes(Mesh * mesh) {
	this->create_node_element_caches();

	ElementDataFieldMap::iterator it=element_data_values_->find("Points");
	ASSERT(it != element_data_values_->end()).error("Missing cache of Points section. Did you call create_node_element_caches()?\n");

	// create nodes of mesh
	std::vector<double> &vect = *( dynamic_cast<ElementDataCache<double> &>(*(it->second)).get_component_data(0).get() );
	unsigned int n_nodes = vect.size()/3;
    mesh->init_node_vector( n_nodes );
	arma::vec3 point;
	for (unsigned int i=0, ivec=0; i<n_nodes; ++i) {
        point(0)=vect[ivec]; ++ivec;
        point(1)=vect[ivec]; ++ivec;
        point(2)=vect[ivec]; ++ivec;
        mesh->add_node(i, point);
	}

	bulk_elements_id_.clear();
}


void VtkMeshReader::read_elements(Mesh * mesh) {
    // read offset section in VTK file
	ElementDataFieldMap::iterator offset_it=element_data_values_->find("offsets");
	ASSERT(offset_it != element_data_values_->end()).error("Missing cache of offsets section. Did you call create_node_element_caches()?\n");
    std::vector<unsigned int> &offsets_vec = *( dynamic_cast<ElementDataCache<unsigned int> &>(*(offset_it->second)).get_component_data(0).get() );

    // read connectivity data section
	ElementDataFieldMap::iterator conn_it=element_data_values_->find("connectivity");
	ASSERT(conn_it != element_data_values_->end()).error("Missing cache of offsets section. Did you call create_node_element_caches()?\n");
    std::vector<unsigned int> &connectivity_vec = *( dynamic_cast<ElementDataCache<unsigned int> &>(*(conn_it->second)).get_component_data(0).get() );

    // iterate trough connectivity data, create bulk elements
    // fill bulk_elements_id_ vector
    bulk_elements_id_.clear();
    bulk_elements_id_.resize(offsets_vec.size());
    mesh->init_element_vector(offsets_vec.size());
    unsigned int i_con = 0, last_offset=0, dim;
    vector<unsigned int> node_list;
    for (unsigned int i=0; i<offsets_vec.size(); ++i) { // iterate trough offset - one value for every element
    	dim = offsets_vec[i] - last_offset - 1; // dimension of vtk element
        for ( ; i_con<offsets_vec[i]; ++i_con ) { // iterate trough all nodes of any element
            node_list.push_back( connectivity_vec[i_con] );
        }
		mesh->add_element(i, dim, dim, 0, node_list); // TODO need to set region_id (3rd parameter, now is created from dim)
        bulk_elements_id_[i] = (LongIdx)i;
        node_list.clear();
        last_offset = offsets_vec[i];
    }

    mesh->create_boundary_elements();
}


void VtkMeshReader::create_node_element_caches() {
	ElementDataFieldMap::iterator it=element_data_values_->find("Points");
	if ( it != element_data_values_->end() ) {
		// prevents repeat call of read_nodes
		return;
	}

	ASSERT( !has_compatible_mesh_ ).error();

	has_compatible_mesh_ = true;

	// read Points data section
	HeaderQuery header_params("Points", 0.0, OutputTime::DiscreteSpace::MESH_DEFINITION);
	auto point_header = this->find_header(header_params);
	bulk_elements_id_.resize(point_header.n_entities);
	for (unsigned int i=0; i<bulk_elements_id_.size(); ++i) bulk_elements_id_[i]=i;
	this->get_element_data<double>(point_header.n_entities, point_header.n_components, false, 0 );

	// read offset data section
	HeaderQuery offsets_params("offsets", 0.0, OutputTime::DiscreteSpace::MESH_DEFINITION);
    auto offset_header = this->find_header(offsets_params);
    for (unsigned int i=bulk_elements_id_.size(); i<offset_header.n_entities; ++i) bulk_elements_id_.push_back(i);
    std::vector<unsigned int> &offsets_vec = *( this->get_element_data<unsigned int>(offset_header.n_entities, offset_header.n_components, false, 0 ) );

    // read connectivity data section
    HeaderQuery con_params("connectivity", 0.0, OutputTime::DiscreteSpace::MESH_DEFINITION);
    auto & con_header = this->find_header(con_params);
    con_header.n_entities = offsets_vec[offsets_vec.size()-1];
    for (unsigned int i=bulk_elements_id_.size(); i<con_header.n_entities; ++i) bulk_elements_id_.push_back(i);
    this->get_element_data<unsigned int>(con_header.n_entities, con_header.n_components, false, 0 );

    has_compatible_mesh_ = false;
	bulk_elements_id_.clear();
}



// explicit instantiation of template methods
template unsigned int read_binary_value<unsigned int>(std::istream &data_stream);
template uint64_t read_binary_value<uint64_t>(std::istream &data_stream);
