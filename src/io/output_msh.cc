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
	.derive_from(OutputFormat::input_type);

// TODO: Remove or adjust the following code.
//	// The variant
//	static Selection variant_sel("GMSH variant");
//	    variant_sel.add_value(OutputMSH::VARIANT_ASCII, "ascii",
//	    		"ASCII variant of GMSH file format");
//	    variant_sel.add_value(OutputMSH::VARIANT_BINARY, "binary",
//	    		"Binary variant of GMSH file format (not supported yet)");
//	    variant_sel.finish();


void OutputMSH::write_msh_header(void)
{
    ofstream &file = this->output_time->get_base_file();

    // Write simple header
    file << "$MeshFormat" << endl;
    file << "2" << " 0 " << sizeof(double) << endl;
    file << "$EndMeshFormat" << endl;
}

void OutputMSH::write_msh_geometry(void)
{
    ofstream &file = this->output_time->get_base_file();
    Mesh* mesh = this->output_time->get_mesh();

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
    ofstream &file = this->output_time->get_base_file();
    Mesh* mesh = this->output_time->get_mesh();
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
    ofstream &file = this->output_time->get_base_file();
    long int item_id = 1;

    /* Set precision to max */
    file.precision(std::numeric_limits<float>::digits10);
    file.precision(std::numeric_limits<double>::digits10);

    /* Write ascii data */
    for(std::vector<boost::any>::iterator item = output_data->data.begin();
            item != output_data->data.end();
            ++item, ++item_id)
    {
        if(item->type() == typeid(double)) {
            file << item_id << " " << boost::any_cast<double>(*item) << std::endl;
        } else if(item->type() == typeid(float)) {
            file << item_id << " " << boost::any_cast<float>(*item) << std::endl;
        } else if(item->type() == typeid(int)) {
            file << item_id << " " << boost::any_cast<int>(*item) << std::endl;
        }
    }
}


void OutputMSH::write_msh_ascii_discont_data(OutputData* output_data)
{
    /* TODO */
}


void OutputMSH::write_msh_node_data(double time, int step)
{
    ofstream &file = this->output_time->get_base_file();
    Mesh *mesh = this->output_time->get_mesh();


    if(this->output_time->node_data.empty() == false) {
        for(vector<OutputData*>::iterator data = this->output_time->node_data.begin();
                    data != this->output_time->node_data.end();
                    ++data)
        {
            file << "$NodeData" << endl;

            file << "1" << endl;     // one string tag
            file << "\"" << (*data)->field->name() << "_[" << (*data)->field->units() <<"]\"" << endl;

            file << "1" << endl;     // one real tag
            file << time << endl;    // first real tag = time

            file << "3" << endl;     // 3 integer tags
            file << step << endl;    // step number (start = 0)
            file << (*data)->field->get_spacedim() << endl;   // number of components
            file << (*data)->data.size() << endl;  // number of values

            this->write_msh_ascii_cont_data(*data);

            file << "$EndNodeData" << endl;
        }
    } else if(this->output_time->corner_data.empty() == false) {
        for(vector<OutputData*>::iterator data = this->output_time->corner_data.begin();
                    data != this->output_time->corner_data.end();
                    ++data)
        {
            file << "$ElementNodeData" << endl;

            file << "1" << endl;     // one string tag
            file << "\"" << (*data)->field->name() << "_[" << (*data)->field->units() <<"]\"" << endl;

            file << "1" << endl;     // one real tag
            file << time << endl;    // first real tag = time

            file << "3" << endl;     // 3 integer tags
            file << step << endl;    // step number (start = 0)
            file << (*data)->field->get_spacedim() << endl;   // number of components
            file << mesh->n_elements() << endl; // number of values

            this->write_msh_ascii_discont_data(*data);

            file << "$EndElementNodeData" << endl;
        }
    }
}

void OutputMSH::write_msh_elem_data(double time, int step)
{
    ofstream &file = this->output_time->get_base_file();

    if(this->output_time->elem_data.empty() == false) {
        for(vector<OutputData*>::iterator data = this->output_time->elem_data.begin();
                    data != this->output_time->elem_data.end();
                    ++data)
        {
            file << "$ElementData" << endl;

            file << "1" << endl;     // one string tag
            file << "\"" << (*data)->field->name() << "_[" << (*data)->field->units() <<"]\"" << endl;

            file << "1" << endl;     // one real tag
            file << time << endl;    // first real tag = time

            file << "3" << endl;     // 3 integer tags
            file << step << endl;    // step number (start = 0)
            file << (*data)->field->get_spacedim() << endl;   // number of components
            file << (*data)->data.size() << endl;  // number of values

            this->write_msh_ascii_cont_data(*data);

            file << "$EndElementData" << endl;
        }
    }
}

int OutputMSH::write_head(void)
{
    xprintf(MsgLog, "%s: Writing output file %s ... ", __func__, this->output_time->get_base_filename().c_str());

    this->write_msh_header();

    this->write_msh_geometry();

    this->write_msh_topology();

    xprintf(MsgLog, "O.K.\n");

    return 1;
}

int OutputMSH::write_data(void)
{
    xprintf(MsgLog, "%s: Writing output file %s ... ", __func__, this->output_time->get_base_filename().c_str());

    // Write header with mesh, when it hasn't been written to output file yet
    if(this->header_written == false) {
        this->write_head();
        this->header_written = true;
    }
        
    this->write_msh_node_data(this->output_time->time, this->output_time->current_step);
    this->write_msh_elem_data(this->output_time->time, this->output_time->current_step);

    // Flush stream to be sure everything is in the file now
    this->output_time->get_base_file().flush();

    xprintf(MsgLog, "O.K.\n");

    return 1;
}

int OutputMSH::write_tail(void)
{
    return 1;
}

OutputMSH::OutputMSH(OutputTime *_output_time)
{
    this->output_time = _output_time;
    this->write_head();
}

OutputMSH::OutputMSH(OutputTime *_output_time, const Input::Record &in_rec)
{
    this->output_time = _output_time;
    this->header_written = false;
}

OutputMSH::~OutputMSH()
{
    this->write_tail();
}


