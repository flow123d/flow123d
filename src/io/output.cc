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
 * @file    output.cc
 * @brief   The functions for all outputs (methods of classes: Output and OutputTime).
 *
 */

#include <string>
#include <petsc.h>

#include "system/xio.h"
#include "io/output.h"
#include "io/output_vtk.h"
#include "io/output_msh.h"
#include "mesh/mesh.h"


OutputData::OutputData(string data_name,
        string data_units,
        int *data_data,
        unsigned int size)
{
    name = new string(data_name); units = new string(data_units);
    data = (void*)data_data;
    type = OUT_ARRAY_INT_SCA;
    comp_num = 1;
    num = size;
}

OutputData::OutputData(string data_name,
        string data_units,
        float *data_data,
        unsigned int size)
{
    name = new string(data_name); units = new string(data_units);
    data = (void*)data_data;
    type = OUT_ARRAY_FLOAT_SCA;
    comp_num = 1;
    num = size;
}

OutputData::OutputData(string data_name,
        string data_units,
        double *data_data,
        unsigned int size)
{
    name = new string(data_name); units = new string(data_units);
    data = (void*)data_data;
    type = OUT_ARRAY_DOUBLE_SCA;
    comp_num = 1;
    num = size;
}

OutputData::OutputData(string data_name,
        string data_units,
        std::vector<int> &data_data)
{
    name = new string(data_name); units = new string(data_units);
    data = (void*)&data_data;
    type = OUT_VECTOR_INT_SCA;
    comp_num = 1;
    num = data_data.size();
}

OutputData::OutputData(string data_name,
        string data_units,
        std::vector< vector<int> > &data_data)
{
    name = new string(data_name); units = new string(data_units);
    data = (void*)&data_data;
    type = OUT_VECTOR_INT_VEC;
    comp_num = 3;
    num = data_data.size();
}

OutputData::OutputData(string data_name,
        string data_units,
        std::vector<float> &data_data)
{
    name = new string(data_name); units = new string(data_units);
    data = (void*)&data_data;
    type = OUT_VECTOR_FLOAT_SCA;
    comp_num = 1;
    num = data_data.size();
}

OutputData::OutputData(string data_name,
        string data_units,
        std::vector< vector<float> > &data_data)
{
    name = new string(data_name); units = new string(data_units);
    data = (void*)&data_data;
    type = OUT_VECTOR_FLOAT_VEC;
    comp_num = 3;
    num = data_data.size();
}

OutputData::OutputData(string data_name,
        string data_units,
        std::vector<double> &data_data)
{
    name = new string(data_name); units = new string(data_units);
    data = (void*)&data_data;
    type = OUT_VECTOR_DOUBLE_SCA;
    comp_num = 1;
    num = data_data.size();
}

OutputData::OutputData(string data_name,
        string data_units,
        std::vector< vector<double> > &data_data)
{
    name = new string(data_name); units = new string(data_units);
    data = (void*)&data_data;
    type = OUT_VECTOR_DOUBLE_VEC;
    comp_num = 3;
    num = data_data.size();
}

OutputData::~OutputData()
{
}

Output::Output(Mesh *_mesh, string fname)
{
    int rank=0;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    /* It's possible now to do output to the file only in the first process */
    if(rank!=0) {
        /* TODO: do something, when support for Parallel VTK is added */
        return;
    }

    base_file = new ofstream;

    base_file->open(fname.c_str());
    INPUT_CHECK( base_file->is_open() , "Can not open output file: %s\n", fname.c_str() );
    xprintf(MsgLog, "Writing flow output file: %s ... \n", fname.c_str());

    base_filename = new string(fname);

    mesh = _mesh;
    node_data = new OutputDataVec;
    elem_data = new OutputDataVec;

    base_filename = new string(fname);

    char *format_name = OptGetStr("Output", "POS_format", "VTK_SERIAL_ASCII");

    if(strcmp(format_name,"ASCII") == 0 || strcmp(format_name,"BIN") == 0) {
        this->file_format = GMSH_MSH_ASCII;
        this->output_msh = new OutputMSH(this);
        this->output_vtk = NULL;
    } else if(strcmp(format_name, "VTK_SERIAL_ASCII") == 0 || strcmp(format_name, "VTK_PARALLEL_ASCII") == 0) {
        this->file_format = VTK_SERIAL_ASCII;
        this->output_msh = NULL;
        this->output_vtk = new OutputVTK(this);
    } else {
        xprintf(Warn,"Unknown output file format: %s.\n", format_name );
        this->file_format = NONE;
        this->output_msh = NULL;
        this->output_vtk = NULL;
    }
}

Output::~Output()
{
    int rank=0;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    /* It's possible now to do output to the file only in the first process */
    if(rank!=0) {
        /* TODO: do something, when support for Parallel VTK is added */
        return;
    }

    if(this->output_msh != NULL) {
        delete this->output_msh;
    }

    if(this->output_vtk != NULL) {
        delete this->output_vtk;
    }

    // Free all reference on node and element data
    if(node_data != NULL) {
        delete node_data;
    }

    if(elem_data != NULL) {
        delete elem_data;
    }

    if(base_filename != NULL) {
        delete base_filename;
    }

    if(base_file != NULL) {
        base_file->close();
        delete base_file;
    }

    xprintf(MsgLog, "O.K.\n");
}

int Output::write_head(void)
{
    switch(this->file_format) {
    case GMSH_MSH_ASCII:
        return this->output_msh->write_head();
    case VTK_SERIAL_ASCII:
        return this->output_vtk->write_head();
    default:
        return 0;
    }
    return 0;
}

int Output::write_tail(void)
{
    switch(this->file_format) {
    case GMSH_MSH_ASCII:
        return this->output_msh->write_tail();
    case VTK_SERIAL_ASCII:
        return this->output_vtk->write_tail();
    default:
        return 0;
    }
    return 0;
}

int Output::write_data()
{
    switch(this->file_format) {
    case GMSH_MSH_ASCII:
        return this->output_msh->write_data();
    case VTK_SERIAL_ASCII:
        return this->output_vtk->write_data();
    default:
        return 0;
    }

    return 0;
}

OutputTime::OutputTime(Mesh *_mesh, string fname)
{
    std::vector<OutputData> *node_data;
    std::vector<OutputData> *elem_data;
    Mesh *mesh = _mesh;
    ofstream *base_file;
    string *base_filename;

    base_file = new ofstream;

    base_file->open(fname.c_str());
    INPUT_CHECK( base_file->is_open() , "Can not open output file: %s\n", fname.c_str() );
    xprintf(MsgLog, "Writing flow output file: %s ... \n", fname.c_str());

    base_filename = new string(fname);

    current_step = 0;

    node_data = new OutputDataVec;
    elem_data = new OutputDataVec;

    set_base_file(base_file);
    set_base_filename(base_filename);
    set_mesh(mesh);
    set_node_data(node_data);
    set_elem_data(elem_data);

    char *format_name = OptGetStr("Output", "POS_format", "VTK_SERIAL_ASCII");

    if(strcmp(format_name,"ASCII") == 0 || strcmp(format_name,"BIN") == 0) {
        this->file_format = GMSH_MSH_ASCII;
        this->output_msh = new OutputMSH(this);
        this->output_vtk = NULL;
    } else if(strcmp(format_name, "VTK_SERIAL_ASCII") == 0 || strcmp(format_name, "VTK_PARALLEL_ASCII") == 0) {
        this->file_format = VTK_SERIAL_ASCII;
        this->output_msh = NULL;
        this->output_vtk = new OutputVTK(this);
    } else {
        xprintf(Warn,"Unknown output file format: %s.\n", format_name );
        this->file_format = NONE;
        this->output_msh = NULL;
        this->output_vtk = NULL;
    }

}

OutputTime::~OutputTime(void)
{
}

int OutputTime::write_data(double time)
{
    int ret = 0;

    switch(this->file_format) {
    case GMSH_MSH_ASCII:
        ret = this->output_msh->write_data(time);
        this->current_step++;
        break;
    case VTK_SERIAL_ASCII:
        ret = this->output_vtk->write_data(time);
        this->current_step++;
        break;
    default:
        break;
    }

    return ret;
}
