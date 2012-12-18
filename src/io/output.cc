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
#include "input/accessors.hh"


using namespace Input::Type;

Record OutputTime::input_type
	= Record("OutputStrem", "Parameters of output.")
	// The name
	.declare_key("name", String(), Default::obligatory(),
			"The name of this stream. Used to reference the output stream.")
	// The stream
	.declare_key("file", FileName::output(), Default::obligatory(),
			"File path to the output stream.")
	// The format
	.declare_key("format", OutputFormat::input_type, Default::optional(),
			"Format of output stream and possible parameters.");


AbstractRecord OutputFormat::input_type
	= AbstractRecord("OutputFormat",
    		"Format of output stream and possible parameters.");
    // Complete declaration of  abstract record OutputFormat




OutputData::OutputData(std::string data_name,
        std::string data_units,
        int *data_data,
        unsigned int size)
{
    name = new string(data_name); units = new string(data_units);
    data = (void*)data_data;
    type = OUT_ARRAY_INT_SCA;
    comp_num = 1;
    num = size;
}

OutputData::OutputData(std::string data_name,
        std::string data_units,
        float *data_data,
        unsigned int size)
{
    name = new string(data_name); units = new string(data_units);
    data = (void*)data_data;
    type = OUT_ARRAY_FLOAT_SCA;
    comp_num = 1;
    num = size;
}

OutputData::OutputData(std::string data_name,
        std::string data_units,
        double *data_data,
        unsigned int size)
{
    name = new string(data_name); units = new string(data_units);
    data = (void*)data_data;
    type = OUT_ARRAY_DOUBLE_SCA;
    comp_num = 1;
    num = size;
}

OutputData::OutputData(std::string data_name,
        std::string data_units,
        std::vector<int> &data_data)
{
    name = new string(data_name); units = new string(data_units);
    data = (void*)&data_data;
    type = OUT_VECTOR_INT_SCA;
    comp_num = 1;
    num = data_data.size();
}

OutputData::OutputData(std::string data_name,
        std::string data_units,
        std::vector< vector<int> > &data_data)
{
    name = new string(data_name); units = new string(data_units);
    data = (void*)&data_data;
    type = OUT_VECTOR_INT_VEC;
    comp_num = 3;
    num = data_data.size();
}

OutputData::OutputData(std::string data_name,
        std::string data_units,
        std::vector<float> &data_data)
{
    name = new string(data_name); units = new string(data_units);
    data = (void*)&data_data;
    type = OUT_VECTOR_FLOAT_SCA;
    comp_num = 1;
    num = data_data.size();
}

OutputData::OutputData(std::string data_name,
        std::string data_units,
        std::vector< vector<float> > &data_data)
{
    name = new string(data_name); units = new string(data_units);
    data = (void*)&data_data;
    type = OUT_VECTOR_FLOAT_VEC;
    comp_num = 3;
    num = data_data.size();
}

OutputData::OutputData(std::string data_name,
        std::string data_units,
        std::vector<double> &data_data)
{
    name = new string(data_name); units = new string(data_units);
    data = (void*)&data_data;
    type = OUT_VECTOR_DOUBLE_SCA;
    comp_num = 1;
    num = data_data.size();
}

OutputData::OutputData(std::string data_name,
        std::string data_units,
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

    // Get number of corners
    unsigned int li;
    this->corner_count = 0;
    FOR_ELEMENTS(mesh, ele) {
        FOR_ELEMENT_NODES(ele, li) {
            this->corner_count++;
        }
    }

    base_filename = new string(fname);

    mesh = _mesh;
    node_data = new OutputDataVec;
    corner_data = new OutputDataVec;
    elem_data = new OutputDataVec;

    base_filename = new string(fname);

    this->file_format = VTK;
    this->output_format = new OutputVTK(this);
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

    if(this->output_format != NULL) {
        delete this->output_format;
    }

    // Free all reference on node and element data
    if(node_data != NULL) {
        delete node_data;
    }

    if(corner_data != NULL) {
        delete corner_data;
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
	if(this->output_format != NULL) {
        return this->output_format->write_head();
    }
    return 0;
}

int Output::write_tail(void)
{
	if(this->output_format != NULL) {
        return this->output_format->write_tail();
    }
    return 0;
}

int Output::write_data()
{
    int rank=0;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if (rank!=0) return 1;

	if(this->output_format != NULL) {
        return this->output_format->write_data();
    }
    return 0;
}

/**
 * \brief This method add right suffix to .pvd VTK file
 */
static inline void fix_VTK_file_name(string *fname)
{
	// When VTK file doesn't .pvd suffix, then add .pvd suffix to this file name
	if(fname->compare(fname->size()-4, 4, ".pvd") != 0) {
		xprintf(Warn, "Renaming name of output file from: %s to %s.pvd\n", fname->c_str(), fname->c_str());
		*fname = *fname + ".pvd";
	}
}

/**
 * \brief This method add right suffix to .msh GMSH file
 */
static inline void fix_GMSH_file_name(string *fname)
{
	// When GMSH file doesn't .msh suffix, then add .msh suffix to this file name
	if(fname->compare(fname->size()-4, 4, ".msh") != 0) {
		xprintf(Warn, "Renaming name of output file from: %s to %s.msh\n", fname->c_str(), fname->c_str());
		*fname = *fname + ".msh";
	}
}



/* Initialize static member of the class */
OutputTime** OutputTime::output_streams = NULL;

/* Initialize static member of the class */
int OutputTime::output_streams_count = 0;


OutputTime *OutputTime::is_created(const Input::Record &in_rec)
{
    string name = in_rec.val<string>("name");

    xprintf(MsgLog, "Trying to find output_stream: %s ... ", name.c_str());

    if(OutputTime::output_streams != NULL) {
        for(int i=0; i<OutputTime::output_streams_count; i++) {
            if(name == *OutputTime::output_streams[i]->name) {
                xprintf(MsgLog, "FOUND\n");
                return OutputTime::output_streams[i];
            }
        }
    }

    xprintf(MsgLog, "NOT FOUND\n");

    return NULL;
}

OutputTime::OutputTime(Mesh *_mesh, const Input::Record &in_rec)
{
    output_format=NULL;
    int rank=0;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if (rank!=0) return;

    OutputTime::output_streams_count++;

    std::vector<OutputData> *node_data;
    std::vector<OutputData> *corner_data;
    std::vector<OutputData> *elem_data;

    Mesh *mesh = _mesh;
    ofstream *base_file;
    string *base_filename;

    string fname = in_rec.val<FilePath>("file");
    string stream_name = in_rec.val<string>("name");

    Input::Iterator<Input::AbstractRecord> format = Input::Record(in_rec).find<Input::AbstractRecord>("format");
    

    // Check if file suffix is suffix of specified file format
    if(format) {
    	if((*format).type() == OutputVTK::input_type) {
    		// This should be pvd file format
    		fix_VTK_file_name(&fname);
    	} else if((*format).type() == OutputMSH::input_type) {
    		// This should be msh file format
    		fix_GMSH_file_name(&fname);
    	} else {
    		// Unsuported file format
    		fix_VTK_file_name(&fname);
    	}
    } else {
		// Default file format is VTK
    	fix_VTK_file_name(&fname);
    }

    base_file = new ofstream;

    base_file->open(fname.c_str());
    INPUT_CHECK( base_file->is_open() , "Can not open output file: %s\n", fname.c_str() );
    xprintf(MsgLog, "Writing flow output file: %s ... \n", fname.c_str());

    // Get number of corners
    unsigned int li, count = 0;
    FOR_ELEMENTS(mesh, ele) {
        FOR_ELEMENT_NODES(ele, li) {
            count++;
        }
    }
    this->set_corner_count(count);

    base_filename = new string(fname);

    this->name = new string(stream_name);
    this->current_step = 0;

    node_data = new OutputDataVec;
    corner_data = new OutputDataVec;
    elem_data = new OutputDataVec;

    set_base_file(base_file);
    set_base_filename(base_filename);
    set_mesh(mesh);
    set_node_data(node_data);
    set_corner_data(corner_data);
    set_elem_data(elem_data);

    if(format) {
		if((*format).type() == OutputVTK::input_type) {
			this->output_format = new OutputVTK(this, *format);
		} else if ( (*format).type() == OutputMSH::input_type) {
			this->output_format = new OutputMSH(this, *format);
		} else {
			xprintf(Warn, "Unsupported file format, using default VTK\n");
			this->output_format = new OutputVTK(this);
		}
    } else {
    	this->output_format = new OutputVTK(this);
    }

}

OutputTime::~OutputTime(void)
{
}

int OutputTime::write_data(double time)
{
    int ret = 0;

    DBGMSG("write_data before MPI rank test\n");
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if (rank != 0 ) return 0;
    DBGMSG("write_data AFTER MPI rank test\n");

    if(this->output_format != NULL) {
    	ret = this->output_format->write_data(time);
    	this->current_step++;
    }

    return ret;
}



