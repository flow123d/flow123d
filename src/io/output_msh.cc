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

void OutputMSH::write_msh_header(void)
{
    ofstream &file = this->output->get_base_file();

    // Write simple header
    file << "$MeshFormat" << endl;
    file << "2" << " 0 " << sizeof(double) << endl;
    file << "$EndMeshFormat" << endl;
}

void OutputMSH::write_msh_geometry(void)
{
    ofstream &file = this->output->get_base_file();
    Mesh* mesh = this->output->get_mesh();

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
    ofstream &file = this->output->get_base_file();
    Mesh* mesh = this->output->get_mesh();
    int i;

    // Write information about elements
    file << "$Elements" << endl;
    file << mesh->n_elements() << endl;
    FOR_ELEMENTS(mesh, elm) {
        // element_id element_type 3_other_tags material region partition
        file << ELEM_FULL_ITER(mesh, elm).index() + 1 << " " << elm->type << " 3 " << elm->mid << " " << elm->rid << " " << elm->pid;

        FOR_ELEMENT_NODES(elm, i)
            file << " " << NODE_FULL_ITER(mesh, elm->node[i]).index() + 1;
        file << endl;
    }
    file << "$EndElements" << endl;
}

void OutputMSH::write_msh_ascii_data(OutputData *out_data)
{
    ofstream &file = this->output->get_base_file();
    long int id = 1;

    switch(out_data->type) {
    case OutputData::OUT_VECTOR_INT_SCA:
        for( std::vector<int>::iterator item = ((std::vector<int>*)out_data->data)->begin();
                item != ((std::vector<int>*)out_data->data)->end();
                ++item, ++id) {
            file << id << " " << *item << endl;
        }
        break;
    case OutputData::OUT_VECTOR_INT_VEC:
        for( std::vector< vector<int> >::iterator vec = ((std::vector< vector<int> >*)out_data->data)->begin();
                vec != ((std::vector< vector<int> >*)out_data->data)->end();
                ++vec, ++id)
        {
            file << id << " ";
            for (std::vector<int>::iterator item = vec->begin();
                    item != vec->end();
                    ++item) {
                file << *item << " ";
            }
            file << endl;
        }
        break;
    case OutputData::OUT_VECTOR_FLOAT_SCA:
        for( std::vector<float>::iterator item = ((std::vector<float>*)out_data->data)->begin();
                item != ((std::vector<float>*)out_data->data)->end();
                ++item, ++id) {
            file << id << " " << *item << endl;
        }
        break;
    case OutputData::OUT_VECTOR_FLOAT_VEC:
        for( std::vector< vector<float> >::iterator vec = ((std::vector< vector<float> >*)out_data->data)->begin();
                vec != ((std::vector< vector<float> >*)out_data->data)->end();
                ++vec)
        {
            file << id << " ";
            for (std::vector<float>::iterator item = vec->begin();
                    item != vec->end();
                    ++item) {
                file << *item << " ";
            }
            file << endl;
        }
        break;
    case OutputData::OUT_VECTOR_DOUBLE_SCA:
        for( std::vector<double>::iterator item = ((std::vector<double>*)out_data->data)->begin();
                item != ((std::vector<double>*)out_data->data)->end();
                ++item, ++id) {
            file << id << " " << *item << endl;
        }
        break;
    case OutputData::OUT_VECTOR_DOUBLE_VEC:
        for( std::vector< vector<double> >::iterator vec = ((std::vector< vector<double> >*)out_data->data)->begin();
                vec != ((std::vector< vector<double> >*)out_data->data)->end();
                ++vec, ++id)
        {
            file << id << " ";
            for (std::vector<double>::iterator item = vec->begin();
                    item != vec->end();
                    ++item) {
                file << *item << " ";
            }
            file << endl;
        }
        break;
    case OutputData::OUT_ARRAY_INT_SCA:
        for(int i=0; i<out_data->num; i++, id++) {
            file << id << " " << ((int*)out_data->data)[i] << endl;
        }
        break;
    case OutputData::OUT_ARRAY_FLOAT_SCA:
        for(int i=0; i<out_data->num; i++, id++) {
            file << id << " " << ((float*)out_data->data)[i] << endl;
        }
        break;
    case OutputData::OUT_ARRAY_DOUBLE_SCA:
        for(int i=0; i<out_data->num; i++, id++) {
            file << id << " " << ((double*)out_data->data)[i] << endl;
        }
        break;
    default:
        xprintf(Err, "This type of data: %d is not supported by MSH file format\n", out_data->type);
        break;
    }
}

void OutputMSH::write_msh_node_data(double time, int step)
{
    ofstream &file = this->output->get_base_file();
    std::vector<OutputData> *node_data = this->output->get_node_data();

    if(node_data != NULL) {
        for(OutputDataVec::iterator dta = node_data->begin();
                    dta != node_data->end();
                    ++dta)
        {
            file << "$NodeData" << endl;

            file << "1" << endl;     // one string tag
            file << "\"" << *dta->getName() << "_[" << *dta->getUnits() <<"]\"" << endl;

            file << "1" << endl;     // one real tag
            file << time << endl;    // first real tag = time

            file << "3" << endl;     // 3 integer tags
            file << step << endl;    // step number (start = 0)
            file << dta->getCompNum() << endl;   // number of components
            file << dta->getValueNum() << endl;  // number of values

            this->write_msh_ascii_data(&(*dta));

            file << "$EndNodeData" << endl;
        }
    }
}

void OutputMSH::write_msh_elem_data(double time, int step)
{
    ofstream &file = this->output->get_base_file();
    std::vector<OutputData> *elem_data = this->output->get_elem_data();

    if(elem_data != NULL) {
        for(OutputDataVec::iterator dta = elem_data->begin();
                    dta != elem_data->end();
                    ++dta)
        {
            file << "$ElementData" << endl;

            file << "1" << endl;     // one string tag
            file << "\"" << *dta->getName() << "_[" << *dta->getUnits() <<"]\"" << endl;

            file << "1" << endl;     // one real tag
            file << time << endl;    // first real tag = time

            file << "3" << endl;     // 3 integer tags
            file << step << endl;    // step number (start = 0)
            file << dta->getCompNum() << endl;   // number of components
            file << dta->getValueNum() << endl;  // number of values

            this->write_msh_ascii_data(&(*dta));

            file << "$EndElementData" << endl;
        }
    }
}

int OutputMSH::write_data(void)
{
    xprintf(MsgLog, "%s: Writing output file %s ... ", __func__, this->output->get_base_filename().c_str());

    this->write_msh_header();

    this->write_msh_geometry();

    this->write_msh_topology();

    this->write_msh_node_data(0.0, 0);

    this->write_msh_elem_data(0.0, 0);

    xprintf(MsgLog, "O.K.\n");

    return 1;
}

int OutputMSH::write_head(void)
{
    xprintf(MsgLog, "%s: Writing output file %s ... ", __func__, this->output->get_base_filename().c_str());

    this->write_msh_header();

    this->write_msh_geometry();

    this->write_msh_topology();

    xprintf(MsgLog, "O.K.\n");

    return 1;
}

int OutputMSH::write_data(double time)
{
    xprintf(MsgLog, "%s: Writing output file %s ... ", __func__, this->output->get_base_filename().c_str());

    this->write_msh_node_data(time, this->output_time->current_step);

    this->write_msh_elem_data(time, this->output_time->current_step);

    // Flush stream to be sure everything is in the file now
    this->output->get_base_file().flush();

    xprintf(MsgLog, "O.K.\n");

    return 1;
}

int OutputMSH::write_tail(void)
{
    return 1;
}

OutputMSH::OutputMSH(Output *_output)
{
    this->output = _output;
    this->write_head();
}

OutputMSH::OutputMSH(OutputTime *_output_time)
{
    this->output = _output_time;
    this->output_time = _output_time;
    this->write_head();
}

OutputMSH::OutputMSH(OutputTime *_output_time, const Input::Record &in_rec)
{
    this->output = _output_time;
    this->output_time = _output_time;
    this->write_head();
}

OutputMSH::~OutputMSH()
{
    this->write_tail();
}

Input::Type::Record & OutputMSH::get_input_type()
{
	using namespace Input::Type;
	static Record rec("OutputMSH", "Parameters of gmsh output format.");

	if (!rec.is_finished()) {

		//rec.derive_from(OutputTime::get_input_type_output_format());

		// The variant
		static Selection variant_sel("GMSH variant");
	    variant_sel.add_value(OutputMSH::VARIANT_ASCII, "ascii",
	    		"ASCII variant of GMSH file format");
	    variant_sel.add_value(OutputMSH::VARIANT_BINARY, "binary",
	    		"Binary variant of GMSH file format (not supported yet)");
	    variant_sel.finish();

		rec.finish();
	}

	return rec;
}
