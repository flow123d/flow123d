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

#include "io/output.h"
#include "io/output_msh.h"
#include "system/xio.h"
#include "mesh/mesh.h"


using namespace Input::Type;

Record OutputMSH::input_type
	= Record("gmsh", "Parameters of gmsh output format.")
	// It is derived from abstract class
	.derive_from(OutputTime::input_format_type);


void OutputMSH::write_msh_header(void)
{
    ofstream &file = this->get_base_file();

    // Write simple header
    file << "$MeshFormat" << endl;
    file << "2" << " 0 " << sizeof(double) << endl;
    file << "$EndMeshFormat" << endl;
}

void OutputMSH::write_msh_geometry(void)
{
    ofstream &file = this->get_base_file();
    Mesh* mesh = this->get_mesh();

    // Write information about nodes
    file << "$Nodes" << endl;
    file <<  mesh->node_vector.size() << endl;
    FOR_NODES(mesh, nod) {
        file << NODE_FULL_ITER(mesh, nod).index() + 1 << " " << nod->getX() << " " << nod->getY() << " " << nod->getZ() << endl;
    }
    file << "$EndNodes" << endl;
}

void OutputMSH::write_msh_topology(void)
{
    ofstream &file = this->get_base_file();
    Mesh* mesh = this->get_mesh();
    unsigned int i;
    const static unsigned int gmsh_simplex_types_[4] = {0, 1, 2, 4};

    // Write information about elements
    file << "$Elements" << endl;
    file << mesh->n_elements() << endl;
    FOR_ELEMENTS(mesh, elm) {
        // element_id element_type 3_other_tags material region partition
        file << ELEM_FULL_ITER(mesh, elm).index() + 1
             << " " << gmsh_simplex_types_[ elm->dim() ]
             << " 3 " << elm->region().id() << " " << elm->region().id() << " " << elm->pid;

        FOR_ELEMENT_NODES(elm, i)
            file << " " << NODE_FULL_ITER(mesh, elm->node[i]).index() + 1;
        file << endl;
    }
    file << "$EndElements" << endl;
}


void OutputMSH::write_msh_ascii_cont_data(OutputData* output_data)
{
    ofstream &file = this->get_base_file();
    long int item_id = 1;
    unsigned int i, j, offset;

    /* Set precision to max */
    file.precision(std::numeric_limits<float>::digits10);
    file.precision(std::numeric_limits<double>::digits10);

    offset = output_data->spacedim;

    /* Write ascii data */
    if(output_data->data_type == OutputData::INT) {
        for(i=0; i < output_data->item_count; i += offset, ++item_id) {
            file << item_id << " ";
            for(j=0; j < offset; j++) {
                file << ((int*)output_data->data)[i*offset + j] << " ";
            }
            file << std::endl;
        }
    } else if(output_data->data_type == OutputData::UINT) {
        for(i=0; i < output_data->item_count; i += offset, ++item_id) {
            file << item_id << " ";
            for(j=0; j < offset; j++) {
                file << ((unsigned int*)output_data->data)[i*offset + j] << " ";
            }
            file << std::endl;
        }
    } else if(output_data->data_type == OutputData::DOUBLE) {
        for(i=0; i < output_data->item_count; i += offset, ++item_id) {
            file << item_id << " ";
            for(j=0; j < offset; j++) {
                file << ((double*)output_data->data)[i*offset + j] << " ";
            }
            file << std::endl;
        }
    }

}


void OutputMSH::write_msh_ascii_discont_data(OutputData* output_data)
{
    /* TODO */
}


void OutputMSH::write_msh_node_data(double time, int step)
{
    ofstream &file = this->get_base_file();
    Mesh *mesh = this->get_mesh();
    OutputData *output_data;

    if(this->node_data.empty() == false) {
        for(vector<OutputData*>::iterator data = this->node_data.begin();
                    data != this->node_data.end();
                    ++data)
        {
            output_data = *data;

            file << "$NodeData" << endl;

            file << "1" << endl;     // one string tag
            file << "\"" << output_data->field->name() << "_[" << output_data->field->units() <<"]\"" << endl;

            file << "1" << endl;     // one real tag
            file << time << endl;    // first real tag = time

            file << "3" << endl;     // 3 integer tags
            file << step << endl;    // step number (start = 0)
            file << output_data->spacedim << endl;   // number of components
            file << output_data->item_count << endl;  // number of values

            this->write_msh_ascii_cont_data(output_data);

            file << "$EndNodeData" << endl;
        }
    } else if(this->corner_data.empty() == false) {
        for(vector<OutputData*>::iterator data = this->corner_data.begin();
                    data != this->corner_data.end();
                    ++data)
        {
            output_data = *data;

            file << "$ElementNodeData" << endl;

            file << "1" << endl;     // one string tag
            file << "\"" << output_data->field->name() << "_[" << output_data->field->units() <<"]\"" << endl;

            file << "1" << endl;     // one real tag
            file << time << endl;    // first real tag = time

            file << "3" << endl;     // 3 integer tags
            file << step << endl;    // step number (start = 0)
            file << output_data->spacedim << endl;   // number of components
            file << mesh->n_elements() << endl; // number of values

            this->write_msh_ascii_discont_data(output_data);

            file << "$EndElementNodeData" << endl;
        }
    }
}

void OutputMSH::write_msh_elem_data(double time, int step)
{
	OutputData* output_data;
    ofstream &file = this->get_base_file();

    if(this->elem_data.empty() == false) {
        for(vector<OutputData*>::iterator data = this->elem_data.begin();
                    data != this->elem_data.end();
                    ++data)
        {
            output_data = *data;
            file << "$ElementData" << endl;

            file << "1" << endl;     // one string tag
            file << "\"" << output_data->field->name() << "_[" << output_data->field->units() <<"]\"" << endl;

            file << "1" << endl;     // one real tag
            file << time << endl;    // first real tag = time

            file << "3" << endl;     // 3 integer tags
            file << step << endl;    // step number (start = 0)
            file << output_data->spacedim << endl;   // number of components
            file << output_data->item_count << endl;  // number of values

            this->write_msh_ascii_cont_data(output_data);

            file << "$EndElementData" << endl;
        }
    }
}

int OutputMSH::write_head(void)
{
    xprintf(MsgLog, "%s: Writing output file %s ... ", __func__, this->base_filename()->c_str());

    this->write_msh_header();

    this->write_msh_geometry();

    this->write_msh_topology();

    xprintf(MsgLog, "O.K.\n");

    return 1;
}

int OutputMSH::write_data(void)
{
    xprintf(MsgLog, "%s: Writing output file %s ... ", __func__, this->base_filename()->c_str());

    // Write header with mesh, when it hasn't been written to output file yet
    if(this->header_written == false) {
        this->write_head();
        this->header_written = true;
    }
        
    this->write_msh_node_data(this->time, this->current_step);
    this->write_msh_elem_data(this->time, this->current_step);

    // Flush stream to be sure everything is in the file now
    this->get_base_file().flush();

    xprintf(MsgLog, "O.K.\n");

    return 1;
}

int OutputMSH::write_tail(void)
{
    return 1;
}

OutputMSH::OutputMSH(const Input::Record &in_rec) : OutputTime(in_rec)
{
	this->file_format = OutputTime::GMSH;
    this->header_written = false;
}

OutputMSH::~OutputMSH()
{
    this->write_tail();
}

