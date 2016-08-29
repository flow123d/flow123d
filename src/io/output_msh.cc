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
#include "mesh/mesh.h"
#include "output_data_base.hh"
#include "input/factory.hh"


FLOW123D_FORCE_LINK_IN_CHILD(gmsh)


using namespace Input::Type;


/**
 * Auxiliary implementation of OutputDataBase that performs output of single zero data for the fields that are
 * off for current time frame.
 */
class DummyOutputData : public OutputDataBase {
public:

    DummyOutputData(std::string field_name_in, OutputDataBase::NumCompValueType n_elem_in)
   {
        this->output_field_name = field_name_in;
        this->n_elem_ = n_elem_in;
        this->n_values = 1;
    }

    virtual ~DummyOutputData() override
    {}

    void print_ascii(ostream &out_stream, unsigned int idx) override
    {
        for(unsigned int i=0; i< n_elem_;i++) out_stream << 0 << " ";
    }

    void print_ascii_all(ostream &out_stream) override
    {
        for(unsigned int i=0; i< n_elem_;i++) out_stream << 0 << " ";
    }

    void print_binary_all(ostream &out_stream) override
    {
        ASSERT(false).error("Not implemented.");
    }

    void print_all_yaml(ostream &out_stream, unsigned int precision) override
    {}

    void get_min_max_range(double &min, double &max) override
    {}
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
    Mesh* mesh = this->_mesh;

    // Write information about nodes
    file << "$Nodes" << endl;
    file <<  mesh->node_vector.size() << endl;
    FOR_NODES(mesh, nod) {
        file << NODE_FULL_ITER(mesh, nod).id() << " " << nod->getX() << " " << nod->getY() << " " << nod->getZ() << endl;
    }
    file << "$EndNodes" << endl;
}

void OutputMSH::write_msh_topology(void)
{
    ofstream &file = this->_base_file;
    Mesh* mesh = this->_mesh;
    unsigned int i;
    const static unsigned int gmsh_simplex_types_[4] = {0, 1, 2, 4};

    // Write information about elements
    file << "$Elements" << endl;
    file << mesh->n_elements() << endl;
    FOR_ELEMENTS(mesh, elm) {
        // element_id element_type 3_other_tags material region partition
        file << ELEM_FULL_ITER(mesh, elm).id()
             << " " << gmsh_simplex_types_[ elm->dim() ]
             << " 3 " << elm->region().id() << " " << elm->region().id() << " " << elm->pid;

        FOR_ELEMENT_NODES(elm, i)
            file << " " << NODE_FULL_ITER(mesh, elm->node[i]).id();
        file << endl;
    }
    file << "$EndElements" << endl;
}


template<class element>
void OutputMSH::write_msh_ascii_cont_data(flow::VectorId<element> &vec, OutputDataPtr output_data)
{
    ofstream &file = this->_base_file;

    /* Set precision to max */
    file.precision(std::numeric_limits<double>::digits10);

    for(unsigned int i=0; i < output_data->n_values; i ++) {
        file << vec(i).id() << " ";
        output_data->print_ascii(file, i);
        file << std::endl;
    }

}


void OutputMSH::write_msh_ascii_discont_data(OutputDataPtr output_data)
{
    Mesh *mesh = this->_mesh;
    ofstream &file = this->_base_file;

    /* Set precision to max */
    file.precision(std::numeric_limits<double>::digits10);

    /* Write ascii data */
    unsigned int i_node;
	unsigned int i_corner = 0;
    FOR_ELEMENTS(mesh, ele) {
        file << ele.id() << " " << ele->n_nodes() << " ";

        FOR_ELEMENT_NODES(ele, i_node) {
            output_data->print_ascii(file, i_corner++);
        }

        file << std::endl;
    }
}



void OutputMSH::write_node_data(OutputDataPtr output_data)
{
    ofstream &file = this->_base_file;
    double time_fixed = isfinite(this->time)?this->time:0;


    file << "$NodeData" << endl;

    file << "1" << endl;     // one string tag
    file << "\"" << output_data->output_field_name <<"\"" << endl;

    file << "1" << endl;     // one real tag
    file << time_fixed << endl;    // first real tag = time

    file << "3" << endl;     // 3 integer tags
    file << this->current_step << endl;    // step number (start = 0)
    file << output_data->n_elem_ << endl;   // number of components
    file << output_data->n_values << endl;  // number of values

    this->write_msh_ascii_cont_data(this->_mesh->node_vector, output_data);

    file << "$EndNodeData" << endl;
}


void OutputMSH::write_corner_data(OutputDataPtr output_data)
{
    ofstream &file = this->_base_file;
    double time_fixed = isfinite(this->time)?this->time:0;

    file << "$ElementNodeData" << endl;

    file << "1" << endl;     // one string tag
    file << "\"" << output_data->output_field_name <<"\"" << endl;

    file << "1" << endl;     // one real tag
    file << time_fixed << endl;    // first real tag = time

    file << "3" << endl;     // 3 integer tags
    file << this->current_step << endl;    // step number (start = 0)
    file << output_data->n_elem_ << endl;   // number of components
    file << this->_mesh->n_elements() << endl; // number of values

    this->write_msh_ascii_discont_data(output_data);

    file << "$EndElementNodeData" << endl;
}

void OutputMSH::write_elem_data(OutputDataPtr output_data)
{
    ofstream &file = this->_base_file;
    double time_fixed = isfinite(this->time)?this->time:0;

    file << "$ElementData" << endl;

    file << "1" << endl;     // one string tag
    file << "\"" << output_data->output_field_name <<"\"" << endl;

    file << "1" << endl;     // one real tag
    file << time_fixed << endl;    // first real tag = time

    file << "3" << endl;     // 3 integer tags
    file << this->current_step << endl;    // step number (start = 0)
    file << output_data->n_elem_ << endl;   // number of components
    file << output_data->n_values << endl;  // number of values

    this->write_msh_ascii_cont_data(this->_mesh->element, output_data);

    file << "$EndElementData" << endl;
}

void OutputMSH::write_field_data(OutputTime::DiscreteSpace type_idx, void (OutputMSH::* format_fce)(OutputDataPtr) )
{
    auto &dummy_data_list = dummy_data_list_[type_idx];
    auto &data_list = this->output_data_vec_[type_idx];

    if (dummy_data_list.size() == 0) {
        // Collect all output fields
        // If more EquationOutput object with different initial times output into same
        // output stream, we may need to possibly update this list on every output frame.
        for(auto out_ptr : data_list)
            dummy_data_list.push_back( std::make_shared<DummyOutputData>(out_ptr->output_field_name, out_ptr->n_elem_));
    }


    auto data_it = data_list.begin();
    for(auto dummy_it = dummy_data_list.begin(); dummy_it != dummy_data_list.end(); ++dummy_it) {
    	DebugOut().fmt("dummy field: {} data field: {}\n", (*dummy_it)->output_field_name, (*data_it)->output_field_name);
        if ((*dummy_it)->output_field_name == (*data_it)->output_field_name) {
            (this->*format_fce)(*data_it); ++data_it;
        } else {
            (this->*format_fce)(*dummy_it);
        }
    }
    ASSERT( data_it ==  data_list.end() )(data_it - data_list.begin())(data_list.size());
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
    // Write header with mesh, when it hasn't been written to output file yet
    if(this->header_written == false) {
        if(this->rank == 0) {
            this->fix_main_file_extension(".msh");
            this->_base_filename.open_stream( this->_base_file );
        }

        this->write_head();
        this->header_written = true;
    }

    LogOut() << __func__ << ": Writing output file " << this->_base_filename << " ... ";


    this->write_field_data(NODE_DATA, &OutputMSH::write_node_data);
    this->write_field_data(CORNER_DATA, &OutputMSH::write_corner_data);
    this->write_field_data(ELEM_DATA, &OutputMSH::write_elem_data);

    // Flush stream to be sure everything is in the file now
    this->_base_file.flush();

    LogOut() << "O.K.";

    return 1;
}



int OutputMSH::write_tail(void)
{
    return 1;
}



