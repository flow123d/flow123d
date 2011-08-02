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

#include "xio.h"
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

void Output::free_data_from_mesh(void)
{
    if(node_scalar != NULL) {
        delete node_scalar->scalars;
        delete node_scalar;
    }

    if(element_scalar != NULL) {
        delete element_scalar->scalars;
        delete element_scalar;
    }

    if(element_vector != NULL) {
        delete element_vector->vectors;
        delete element_vector;
    }
}

void Output::get_data_from_mesh(void)
{
    int rank=0;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    /* It's possible now to do output to the file only in the first process */
    if(rank!=0) {
        /* TODO: do something, when support for Parallel VTK is added */
        return;
    }

    NodeIter node;
    ElementIter ele;

    node_scalar = new OutScalar;
    element_scalar = new OutScalar;
    element_vector = new OutVector;

    /* Fill temporary vector of node scalars */
    node_scalar->scalars = new ScalarFloatVector;
    node_scalar->name = "node_scalars";
    node_scalar->unit = "";
    /* Generate vector for scalar data of nodes */
    node_scalar->scalars->reserve(mesh->node_vector.size());   // reserver memory for vector
    FOR_NODES( node ) {
        node_scalar->scalars->push_back(node->scalar);
    }

    /* Fill vectors of element scalars and vectors */
    element_scalar->scalars = new ScalarFloatVector;
    element_scalar->name = "element_scalars";
    element_scalar->unit = "";
    element_vector->vectors = new VectorFloatVector;
    element_vector->name = "element_vectors";
    element_vector->unit = "";
    /* Generate vectors for scalar and vector data of nodes */
    element_scalar->scalars->reserve(mesh->n_elements());
    element_vector->vectors->reserve(mesh->n_elements());
    FOR_ELEMENTS(ele) {
        /* Add scalar */
        element_scalar->scalars->push_back(ele->scalar);
        /* Add vector */
        vector<double> vec;
        vec.reserve(3);
        vec.push_back(ele->vector[0]);
        vec.push_back(ele->vector[1]);
        vec.push_back(ele->vector[2]);
        element_vector->vectors->push_back(vec);
    }

    register_node_data(node_scalar->name, node_scalar->unit, *node_scalar->scalars);
    register_elem_data(element_scalar->name, element_scalar->unit, *element_scalar->scalars);
    register_elem_data(element_vector->name, element_vector->unit, *element_vector->vectors);
}

/**
 * \brief NULL function for not yet supported formats
 *
 * \param[in]   *output The pointer at output object.
 *
 * \return This function returns always 0.
 */
int write_null_data(Output *output)
{
    xprintf(Msg, "This file format is not yet supported\n");

    return 0;
}

int Output::write_data(void)
{
    int rank=0;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    /* It's possible now to do output to the file only in the first process */
    if(rank!=0) {
        /* TODO: do something, when support for Parallel VTK is added */
        return 0;
    }

    return _write_data(this);
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

    if( OptGetBool("Output", "Write_output_file", "no") == false ) {
        base_filename = NULL;
        base_file = NULL;
        mesh = NULL;

        return;
    }

    base_file = new ofstream;

    base_file->open(fname.c_str());
    if(base_file->is_open() == false) {
        xprintf(Msg, "Could not write output to the file: %s\n", fname.c_str());
        base_filename = NULL;
        delete base_file;
        base_file = NULL;
        mesh = NULL;

        return;
    } else {
        xprintf(Msg, "Writing flow output file: %s ... \n", fname.c_str());
    }

    base_filename = new string(fname);

    mesh = _mesh;
    node_data = new OutputDataVec;
    elem_data = new OutputDataVec;


    format_type =  parse_output_format( OptGetStr("Output", "POS_format", "VTK_SERIAL_ASCII") );

    switch(format_type) {
    case VTK_SERIAL_ASCII:
    case VTK_PARALLEL_ASCII:
        _write_data = write_vtk_vtu_data;
        break;
    case GMSH_MSH_ASCII:
    case GMSH_MSH_BIN:
        _write_data = write_msh_data;
        break;
    default:
        _write_data = write_null_data;
        break;
    }

    base_filename = new string(fname);
}

OutFileFormat Output::parse_output_format(char* format_name)
{
    if(strcmp(format_name,"ASCII") == 0)
        return GMSH_MSH_ASCII;
    if(strcmp(format_name,"BIN") == 0)
        return GMSH_MSH_BIN;
    if(strcmp(format_name, "VTK_SERIAL_ASCII") == 0)
        return VTK_SERIAL_ASCII;
    if(strcmp(format_name, "VTK_PARALLEL_ASCII") == 0)
        return VTK_PARALLEL_ASCII;
    xprintf(Warn,"Unknown output file format: %s.\n", format_name );
    return (VTK_SERIAL_ASCII);
}


Output::~Output()
{
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

    xprintf(Msg, "O.K.\n");
}

/**
 * \brief This is fake output function for not supported formats. It writes
 * only warning to the stdout and log file.
 *
 * \param[in] *output   The pointer at OutputTime function
 *
 * \return This function always return 0.
 */
int write_null_head(OutputTime *output)
{
    xprintf(Msg, "This file format: %d is not yet supported\n", output->get_format_type());

    return 0;
}

/**
 * \brief This is fake output function for not supported formats. It writes
 * only warning to the stdout and log file.
 *
 * \param[in] *output   The pointer at OutputTime function
 *
 * \return This function always return 0.
 */
int write_null_time_data(OutputTime *output, double time, int step)
{
    xprintf(Msg, "This file format: %d is not yet supported\n", output->get_format_type());

    return 0;
}

/**
 * \brief This is fake output function for not supported formats. It writes
 * only warning to the stdout and log file.
 *
 * \param[in] *output   The pointer at OutputTime function
 *
 * \return This function always return 0.
 */
int write_null_tail(OutputTime *output)
{
    xprintf(Msg, "This file format: %d is not yet supported\n", output->get_format_type());

    return 0;
}

int OutputTime::write_data(double time)
{
    int rank=0;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    /* It's possible now to do output to the file only in the first process */
    if(rank!=0) {
        /* TODO: do something, when support for Parallel VTK is added */
        return 0;
    }
    return _write_data(this, time, current_step++);
}

OutputTime::OutputTime(Mesh *_mesh, string fname)
{
    int rank=0;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    /* It's possible now to do output to the file only in the first process */
    if(rank!=0) {
        /* TODO: do something, when support for Parallel VTK is added */
        return;
    }

    std::vector<OutputData> *node_data;
    std::vector<OutputData> *elem_data;
    Mesh *mesh = _mesh;
    ofstream *base_file;
    string *base_filename;
    OutFileFormat format_type;

    if( OptGetBool("Output", "Write_output_file", "no") == false ) {
        base_filename = NULL;
        base_file = NULL;
        mesh = NULL;

        return;
    }

    base_file = new ofstream;

    base_file->open(fname.c_str());
    if(base_file->is_open() == false) {
        xprintf(Msg, "Could not write output to the file: %s\n", fname.c_str());
        base_filename = NULL;
        delete base_file;
        base_file = NULL;
        mesh = NULL;

        return;
    } else {
        xprintf(Msg, "Writing flow output file: %s ... \n", fname.c_str());
    }

    base_filename = new string(fname);

    current_step = 0;

    node_data = new OutputDataVec;
    elem_data = new OutputDataVec;

    // TODO: remove in the future
    format_type =  parse_output_format( OptGetStr("Output", "POS_format", "VTK_SERIAL_ASCII") );

    set_base_file(base_file);
    set_base_filename(base_filename);
    set_mesh(mesh);
    set_node_data(node_data);
    set_elem_data(elem_data);

    set_format_type(format_type);

    switch(format_type) {
    case VTK_SERIAL_ASCII:
    case VTK_PARALLEL_ASCII:
        _write_head = write_vtk_pvd_head;
        _write_data = write_vtk_pvd_data;
        _write_tail = write_vtk_pvd_tail;
        break;
    case GMSH_MSH_ASCII:
    case GMSH_MSH_BIN:
        _write_head = write_msh_head;
        _write_data = write_msh_time_data;
        _write_tail = write_msh_tail;
        break;
    default:
        _write_head = write_null_head;
        _write_data = write_null_time_data;
        _write_tail = write_null_tail;
        break;
    }

    _write_head(this);
}

OutputTime::~OutputTime(void)
{
    _write_tail(this);
}
