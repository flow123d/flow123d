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
 * $Id: output_msh.cc 2505 2013-09-13 14:52:27Z jiri.hnidek $
 * $Revision$
 * $LastChangedBy$
 * $LastChangedDate$
 *
 * @file    output_msh.cc
 * @brief   The functions for outputs to GMSH files.
 *
 */

#include "output_msh.hh"
#include "mesh/mesh.h"
#include "output_data_base.hh"
#include "input/factory.hh"


FLOW123D_FORCE_LINK_IN_CHILD(gmsh)


using namespace Input::Type;

const Record & OutputMSH::get_input_type() {
	return Record("gmsh", "Parameters of gmsh output format.")
		// It is derived from abstract class
		.derive_from(OutputTime::get_input_format_type())
		.close();
}

const int OutputMSH::registrar = Input::register_class< OutputMSH, const Input::Record & >("gmsh") +
		OutputMSH::get_input_type().size();


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
        output_data->print(file, i);
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
            output_data->print(file, i_corner++);
        }

        file << std::endl;
    }
}


void OutputMSH::write_msh_node_data(double time, int step)
{
    ofstream &file = this->_base_file;
    Mesh *mesh = this->_mesh;

    double time_fixed = isfinite(time)?time:0;

    for(OutputDataPtr output_data :  this->output_data_vec_[NODE_DATA])
        {
            file << "$NodeData" << endl;

            file << "1" << endl;     // one string tag
            file << "\"" << output_data->output_field_name <<"\"" << endl;

            file << "1" << endl;     // one real tag
            file << time_fixed << endl;    // first real tag = time

            file << "3" << endl;     // 3 integer tags
            file << step << endl;    // step number (start = 0)
            file << output_data->n_elem_ << endl;   // number of components
            file << output_data->n_values << endl;  // number of values

            this->write_msh_ascii_cont_data(mesh->node_vector, output_data);

            file << "$EndNodeData" << endl;
        }
    for(OutputDataPtr output_data :  this->output_data_vec_[CORNER_DATA] )
        {
            file << "$ElementNodeData" << endl;

            file << "1" << endl;     // one string tag
            file << "\"" << output_data->output_field_name <<"\"" << endl;

            file << "1" << endl;     // one real tag
            file << time_fixed << endl;    // first real tag = time

            file << "3" << endl;     // 3 integer tags
            file << step << endl;    // step number (start = 0)
            file << output_data->n_elem_ << endl;   // number of components
            file << mesh->n_elements() << endl; // number of values

            this->write_msh_ascii_discont_data(output_data);

            file << "$EndElementNodeData" << endl;
        }
}

void OutputMSH::write_msh_elem_data(double time, int step)
{
    ofstream &file = this->_base_file;

    double time_fixed = isfinite(time) ? time : 0;

    for(OutputDataPtr output_data :  this->output_data_vec_[ELEM_DATA] )
        {
            file << "$ElementData" << endl;

            file << "1" << endl;     // one string tag
            file << "\"" << output_data->output_field_name <<"\"" << endl;

            file << "1" << endl;     // one real tag
            file << time_fixed << endl;    // first real tag = time

            file << "3" << endl;     // 3 integer tags
            file << step << endl;    // step number (start = 0)
            file << output_data->n_elem_ << endl;   // number of components
            file << output_data->n_values << endl;  // number of values

            this->write_msh_ascii_cont_data(this->_mesh->element, output_data);

            file << "$EndElementData" << endl;
        }
}

int OutputMSH::write_head(void)
{
    xprintf(MsgLog, "%s: Writing output file %s ... ", __func__,
            this->_base_filename.c_str());

    this->write_msh_header();

    this->write_msh_geometry();

    this->write_msh_topology();

    xprintf(MsgLog, "O.K.\n");

    return 1;
}

int OutputMSH::write_data(void)
{
    xprintf(MsgLog, "%s: Writing output file %s ... ", __func__,
            this->_base_filename.c_str());

    // Write header with mesh, when it hasn't been written to output file yet
    if(this->header_written == false) {
        this->write_head();
        this->header_written = true;
    }

    this->write_msh_node_data(this->time, this->current_step);
    this->write_msh_elem_data(this->time, this->current_step);

    // Flush stream to be sure everything is in the file now
    this->_base_file.flush();

    xprintf(MsgLog, "O.K.\n");

    return 1;
}



int OutputMSH::write_tail(void)
{
    return 1;
}



OutputMSH::OutputMSH(const Input::Record &in_rec) : OutputTime(in_rec)
{
	this->fix_main_file_extension(".msh");
    this->header_written = false;

    if(this->rank == 0) {
        this->_base_file.open(this->_base_filename.c_str());
        INPUT_CHECK( this->_base_file.is_open() , "Can not open output file: %s\n", this->_base_filename.c_str() );
        xprintf(MsgLog, "Writing flow output file: %s ... \n", this->_base_filename.c_str());
    }
}

OutputMSH::~OutputMSH()
{
    this->write_tail();
}

