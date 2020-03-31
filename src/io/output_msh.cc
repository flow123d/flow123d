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
 * @file    output_msh.cc
 * @brief   The functions for outputs to GMSH files.
 */

#include "output_msh.hh"
#include "output_mesh.hh"
#include "output_element.hh"
#include "mesh/mesh.h"
#include "element_data_cache_base.hh"
#include "input/factory.hh"
#include "tools/unit_si.hh"


FLOW123D_FORCE_LINK_IN_CHILD(gmsh)


using namespace Input::Type;


/**
 * Auxiliary implementation of ElementDataCacheBase that performs output of single zero data for the fields that are
 * off for current time frame.
 */
class DummyOutputData : public ElementDataCacheBase {
public:

    DummyOutputData(std::string field_name_in, unsigned int n_comp_in)
   {
        this->field_input_name_ = field_name_in;
        this->n_comp_ = n_comp_in;
        this->n_values_ = 1;
    }

    virtual ~DummyOutputData() override
    {}

    void print_ascii(ostream &out_stream, unsigned int) override
    {
        for(unsigned int i=0; i< n_comp_;i++) out_stream << 0 << " ";
    }

    void print_ascii_all(ostream &out_stream) override
    {
        for(unsigned int i=0; i< n_comp_;i++) out_stream << 0 << " ";
    }

    void print_binary_all(ostream &, bool) override
    {
        ASSERT(false).error("Not implemented.");
    }

    void print_yaml_subarray(ostream &, unsigned int, unsigned int , unsigned int) override
    {}

    void get_min_max_range(double &, double &) override
    {}

    void read_ascii_data(Tokenizer &, unsigned int, unsigned int ) override
    {}

    void read_binary_data(std::istream &, unsigned int, unsigned int) override
    {}

    std::shared_ptr< ElementDataCacheBase > gather(Distribution *, LongIdx *) override
    {
    	return std::make_shared<DummyOutputData>(this->field_input_name_, this->n_comp_);
    }

    std::shared_ptr< ElementDataCacheBase > element_node_cache_fixed_size(std::vector<unsigned int> &) override
    {
    	return std::make_shared<DummyOutputData>(this->field_input_name_, this->n_comp_);
    }

    std::shared_ptr< ElementDataCacheBase > element_node_cache_optimize_size(std::vector<unsigned int> &) override
    {
    	return std::make_shared<DummyOutputData>(this->field_input_name_, this->n_comp_);
    }

    std::shared_ptr< ElementDataCacheBase > compute_node_data(std::vector<unsigned int> &, unsigned int ) override
    {
    	return std::make_shared<DummyOutputData>(this->field_input_name_, this->n_comp_);
    }

};







const Record & OutputMSH::get_input_type() {
	return Record("gmsh", "Parameters of gmsh output format.")
		// It is derived from abstract class
		.derive_from(OutputTime::get_input_format_type())
		.close();
}

const int OutputMSH::registrar = Input::register_class< OutputMSH >("gmsh") +
		OutputMSH::get_input_type().size();


OutputMSH::OutputMSH()
{
    this->enable_refinement_ = false;
    this->header_written = false;

    dummy_data_list_.resize(OutputTime::N_DISCRETE_SPACES);


}

OutputMSH::~OutputMSH()
{
    this->write_tail();
}




void OutputMSH::write_msh_header(void)
{
    ofstream &file = this->_base_file;

    // Write simple header
    file << "$MeshFormat" << endl;
    file << "2" << " 0 " << sizeof(double) << endl;
    file << "$EndMeshFormat" << endl;
}

void OutputMSH::write_msh_geometry(void)
{
    ofstream &file = this->_base_file;

    // Write information about nodes
    file << "$Nodes" << endl;
    file << this->nodes_->n_values() << endl;
    auto &id_node_vec = *( this->node_ids_->get_component_data(0).get() );
    for(unsigned int i_node=0; i_node < id_node_vec.size(); ++i_node) {
        file << id_node_vec[i_node] << " ";
        this->nodes_->print_ascii(file, i_node);
        file << endl;
    }
    file << "$EndNodes" << endl;
}

void OutputMSH::write_msh_topology(void)
{
    ofstream &file = this->_base_file;
    const static unsigned int gmsh_simplex_types_[4] = {0, 1, 2, 4};
    auto &id_elem_vec = *( this->elem_ids_->get_component_data(0).get() );
    auto &id_node_vec = *( this->node_ids_->get_component_data(0).get() );
    auto &connectivity_vec = *( this->connectivity_->get_component_data(0).get() );
    auto &offsets_vec = *( this->offsets_->get_component_data(0).get() );
    auto &regions_vec = *( this->region_ids_->get_component_data(0).get() );
    auto &partition_vec = *( this->partitions_->get_component_data(0).get() );

    unsigned int n_nodes, i_node=0;

    // Write information about elements
    file << "$Elements" << endl;
    file << this->offsets_->n_values() << endl;
    ElementAccessor<OutputElement::spacedim> elm;
    for(unsigned int i_elm=0; i_elm < id_elem_vec.size(); ++i_elm) {
    	n_nodes = (i_elm==0) ? (offsets_vec[0]) : (offsets_vec[i_elm]-offsets_vec[i_elm-1]);
        // element_id element_type 3_other_tags material region partition
        file << id_elem_vec[i_elm]
             << " " << gmsh_simplex_types_[ n_nodes-1 ]
             << " 3 " << regions_vec[i_elm] << " " << regions_vec[i_elm] << " " << partition_vec[i_elm];

        for(unsigned int i=0; i<n_nodes; i++, i_node++) {
            file << " " << id_node_vec[connectivity_vec[i_node]];
        }
        file << endl;
    }
    file << "$EndElements" << endl;
}


void OutputMSH::write_msh_ascii_data(std::shared_ptr<ElementDataCache<unsigned int>> id_cache, OutputDataPtr output_data, bool discont)
{
	ofstream &file = this->_base_file;
    auto &id_vec = *( id_cache->get_component_data(0).get() );

    if (discont) { // corner data
    	auto &offsets_vec = *( this->offsets_->get_component_data(0).get() );
    	unsigned int n_nodes, i_corner=0;
        for(unsigned int i=0; i < id_vec.size(); ++i) {
        	n_nodes = (i==0) ? (offsets_vec[0]) : (offsets_vec[i]-offsets_vec[i-1]);
            file << id_vec[i] << " " << n_nodes << " ";
            for (unsigned int j=0; j<n_nodes; j++)
            	output_data->print_ascii(file, i_corner++);
            file << std::endl;
        }
    } else { // element / node data
        for(unsigned int i=0; i < output_data->n_values(); ++i) {
            file << id_vec[i] << " ";
            output_data->print_ascii(file, i);
            file << std::endl;
        }

    }
}


void OutputMSH::write_node_data(OutputDataPtr output_data)
{
    ofstream &file = this->_base_file;
    double time_fixed = isfinite(this->time)?this->time:0;
    time_fixed /= UnitSI().s().convert_unit_from(this->unit_string_);

    file << "$NodeData" << endl;

    file << "1" << endl;     // one string tag
    file << "\"" << output_data->field_input_name() <<"\"" << endl;

    file << "1" << endl;     // one real tag
    file << time_fixed << endl;    // first real tag = time

    file << "3" << endl;     // 3 integer tags
    file << this->current_step << endl;    // step number (start = 0)
    file << output_data->n_comp() << endl;   // number of components
    file << output_data->n_values() << endl;  // number of values

    this->write_msh_ascii_data(this->node_ids_, output_data);

    file << "$EndNodeData" << endl;
}


void OutputMSH::write_corner_data(OutputDataPtr output_data)
{
    ofstream &file = this->_base_file;
    double time_fixed = isfinite(this->time)?this->time:0;

    file << "$ElementNodeData" << endl;

    file << "1" << endl;     // one string tag
    file << "\"" << output_data->field_input_name() <<"\"" << endl;

    file << "1" << endl;     // one real tag
    file << time_fixed << endl;    // first real tag = time

    file << "3" << endl;     // 3 integer tags
    file << this->current_step << endl;    // step number (start = 0)
    file << output_data->n_comp() << endl;   // number of components
    file << this->offsets_->n_values() << endl; // number of values

    this->write_msh_ascii_data(this->elem_ids_, output_data, true);

    file << "$EndElementNodeData" << endl;
}

void OutputMSH::write_elem_data(OutputDataPtr output_data)
{
    ofstream &file = this->_base_file;
    double time_fixed = isfinite(this->time)?this->time:0;

    file << "$ElementData" << endl;

    file << "1" << endl;     // one string tag
    file << "\"" << output_data->field_input_name() <<"\"" << endl;

    file << "1" << endl;     // one real tag
    file << time_fixed << endl;    // first real tag = time

    file << "3" << endl;     // 3 integer tags
    file << this->current_step << endl;    // step number (start = 0)
    file << output_data->n_comp() << endl;   // number of components
    file << output_data->n_values() << endl;  // number of values

    this->write_msh_ascii_data(this->elem_ids_, output_data);

    file << "$EndElementData" << endl;
}

int OutputMSH::write_head(void)
{
	LogOut() << __func__ << ": Writing output file " << this->_base_filename << " ... ";

    this->write_msh_header();

    this->write_msh_geometry();

    this->write_msh_topology();

    LogOut() << "O.K.";

    return 1;
}

int OutputMSH::write_data(void)
{
    /* Output of serial format is implemented only in the first process */
    if (this->rank_ != 0) {
        return 0;
    }

    // Write header with mesh, when it hasn't been written to output file yet
    if(this->header_written == false) {
        this->fix_main_file_extension(".msh");
        try {
            this->_base_filename.open_stream( this->_base_file );
            this->set_stream_precision(this->_base_file);
        } INPUT_CATCH(FilePath::ExcFileOpen, FilePath::EI_Address_String, input_record_)

        this->write_head();
        this->header_written = true;
    }

    LogOut() << __func__ << ": Writing output file " << this->_base_filename << " ... ";


    auto &node_data_list = this->output_data_vec_[NODE_DATA];
    for(auto data_it = node_data_list.begin(); data_it != node_data_list.end(); ++data_it) {
    	write_node_data(*data_it);
    }
    auto &corner_data_list = this->output_data_vec_[CORNER_DATA];
    for(auto data_it = corner_data_list.begin(); data_it != corner_data_list.end(); ++data_it) {
    	write_corner_data(*data_it);
    }
    auto &elem_data_list = this->output_data_vec_[ELEM_DATA];
    for(auto data_it = elem_data_list.begin(); data_it != elem_data_list.end(); ++data_it) {
    	write_elem_data(*data_it);
    }

    // Flush stream to be sure everything is in the file now
    this->_base_file.flush();

    LogOut() << "O.K.";

    return 1;
}



int OutputMSH::write_tail(void)
{
    return 1;
}



void OutputMSH::add_dummy_fields()
{
	const std::vector<OutputTime::DiscreteSpace> space_types = {OutputTime::NODE_DATA, OutputTime::CORNER_DATA, OutputTime::ELEM_DATA};
	for (auto type_idx : space_types) {
	    auto &dummy_data_list = dummy_data_list_[type_idx];
	    auto &data_list = this->output_data_vec_[type_idx];

        // Collect all output fields
		if (dummy_data_list.size() == 0)
			for(auto out_ptr : data_list)
				dummy_data_list.push_back( std::make_shared<DummyOutputData>(out_ptr->field_input_name(), out_ptr->n_comp()));

	    auto data_it = data_list.begin();
	    for(auto dummy_it = dummy_data_list.begin(); dummy_it != dummy_data_list.end(); ++dummy_it) {
	        if ( data_it == data_list.end() ) {
	        	data_list.push_back( *dummy_it );
	        } else if ((*dummy_it)->field_input_name() == (*data_it)->field_input_name()) {
	        	++data_it;
	        } else {
	        	data_list.push_back( *dummy_it );
	        }
	    }
	}
}



void OutputMSH::set_output_data_caches(std::shared_ptr<OutputMeshBase> mesh_ptr) {
    OutputTime::set_output_data_caches(mesh_ptr);

    mesh_ptr->get_master_mesh()->create_id_caches();
    this->node_ids_ = mesh_ptr->get_master_mesh()->node_ids_;
    this->elem_ids_ = mesh_ptr->get_master_mesh()->elem_ids_;
    this->region_ids_ = mesh_ptr->get_master_mesh()->region_ids_;
    this->partitions_ = mesh_ptr->get_master_mesh()->partitions_;
}

