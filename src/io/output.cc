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
 * $Id: output.cc 2505 2013-09-13 14:52:27Z jiri.hnidek $
 * $Revision: 2505 $
 * $LastChangedBy: jiri.hnidek $
 * $LastChangedDate: 2013-09-13 16:52:27 +0200 (Pá, 13 IX 2013) $
 *
 * @file    output.cc
 * @brief   The functions for all outputs (methods of classes: Output and OutputTime).
 *
 */

#include <string>
#include <typeinfo>
#include <petsc.h>
#include <boost/any.hpp>
#include <assert.h>

#include "system/xio.h"
#include "io/output.h"
#include "io/output_vtk.h"
#include "io/output_msh.h"
#include "mesh/mesh.h"
#include "input/accessors.hh"


using namespace Input::Type;

Record OutputTime::input_type
    = Record("OutputStream", "Parameters of output.")
    // The name
    .declare_key("name", String(), Default::obligatory(),
            "The name of this stream. Used to reference the output stream.")
    // The stream
    .declare_key("file", FileName::output(), Default::obligatory(),
            "File path to the connected output file.")
            // The format
	.declare_key("format", OutputTime::input_format_type, Default::optional(),
			"Format of output stream and possible parameters.");
#if 0
    // The format
    .declare_key("format", OutputFormat::input_type, Default::optional(),
            "Format of output stream and possible parameters.");

AbstractRecord OutputFormat::input_type
    = AbstractRecord("OutputFormat",
            "Format of output stream and possible parameters.");
    // Complete declaration of  abstract record OutputFormat
#endif

AbstractRecord OutputTime::input_format_type
    = AbstractRecord("OutputTime",
            "Format of output stream and possible parameters.");

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


OutputDataBase *OutputTime::output_data_by_field_name
		(const std::string &multi_field_name, const std::string &field_name, DiscreteSpace ref_type)
{
    std::vector<OutputDataBase*> *data_vector;

    switch(ref_type) {
    case NODE_DATA:
        data_vector = &this->node_data;
        break;
    case CORNER_DATA:
        data_vector = &this->corner_data;
        break;
    case ELEM_DATA:
        data_vector = &this->elem_data;
        break;
    }

    /* Try to find existing data */
    for(auto &data : *data_vector)
        if ((data->field_name == field_name) && (data->multi_field_name == multi_field_name) )     return data;

    return nullptr;
}

/* Initialize static member of the class */
std::vector<OutputTime*> OutputTime::output_streams;

// Destroy all objects
void OutputTime::destroy_all(void)
{
    // Delete all objects
    for(std::vector<OutputTime*>::iterator ot_iter = OutputTime::output_streams.begin();
        ot_iter != OutputTime::output_streams.end();
        ++ot_iter)
    {
        delete *ot_iter;
    }

    OutputTime::output_streams.clear();
}

OutputTime *OutputTime::output_stream_by_name(string name)
{
	OutputTime *output_time;
    // Try to find existing object
    for(std::vector<OutputTime*>::iterator output_iter = OutputTime::output_streams.begin();
            output_iter != OutputTime::output_streams.end();
            ++output_iter)
    {
        output_time = (*output_iter);
        if( *(output_time->name) == name) {
            return output_time;
        }
    }

    return NULL;
}


OutputTime *OutputTime::output_stream_by_key_name(const Input::Record &in_rec, const string key_name)
{
	// TODO: do not try to find empty string and raise exception

	// Try to find record with output stream (the key is name of data)
	Input::Iterator<string> stream_name_iter = in_rec.find<string>(key_name);

	// If record was not found, then throw exception
	if(!stream_name_iter) {
		return nullptr;
	}

	// Try to find existing output stream
	return output_stream_by_name(*stream_name_iter);
}


OutputTime* OutputTime::create_output_stream(const Input::Record &in_rec)
{
    OutputTime* output_time;

    Input::Iterator<Input::AbstractRecord> format = Input::Record(in_rec).find<Input::AbstractRecord>("format");

    if(format) {
        if((*format).type() == OutputVTK::input_type) {
            output_time = new OutputVTK(in_rec);
        } else if ( (*format).type() == OutputMSH::input_type) {
            output_time = new OutputMSH(in_rec);
        } else {
            xprintf(Warn, "Unsupported file format, using default VTK\n");
            output_time = new OutputVTK(in_rec);
        }
    } else {
        output_time = new OutputVTK(in_rec);
    }

    return output_time;
}


OutputTime *OutputTime::output_stream(const Input::Record &in_rec)
{
    // testing rank of process
    int ierr, rank;
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    ASSERT(ierr == 0, "Error in MPI_Comm_rank.");

    /* It's possible now to do output to the file only in the first process */
    if(rank != 0) {
        xprintf(MsgLog, "NOT MASTER PROC\n");
        /* TODO: do something, when support for Parallel VTK is added */
        return NULL;
    }

    OutputTime *output_time;
    string name = in_rec.val<string>("name");

    xprintf(MsgLog, "Trying to find output_stream: %s ... ", name.c_str());

    output_time = OutputTime::output_stream_by_name(name);
    if(output_time != NULL) {
        xprintf(MsgLog, "FOUND\n");
        return output_time;
    }

    xprintf(MsgLog, "NOT FOUND. Creating new ... ");

    output_time = OutputTime::create_output_stream(in_rec);
    OutputTime::output_streams.push_back(output_time);

    xprintf(MsgLog, "DONE\n");

    return output_time;
}


OutputTime::OutputTime(const Input::Record &in_rec)
{
    int ierr;
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    ASSERT(ierr == 0, "Error in MPI_Comm_rank.");
    
    /* It's possible now to do output to the file only in the first process */
    if(rank!=0) {
        /* TODO: do something, when support for Parallel VTK is added */
        return;
    }

    Mesh *mesh = NULL;  // This is set, when first register_* method is called

    ofstream *base_file;
    string *base_filename;

    string fname = in_rec.val<FilePath>("file");
    string stream_name = in_rec.val<string>("name");

    Input::Iterator<Input::AbstractRecord> format = Input::Record(in_rec).find<Input::AbstractRecord>("format");

    // TODO: move this part to OutputVTK.cc and OutputMSH.cc
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

    base_filename = new string(fname);

    this->name = new string(stream_name);
    this->current_step = 0;

    set_base_file(base_file);
    this->_base_filename = base_filename;
    set_mesh(mesh);

    this->time = -1.0;
    this->write_time = -1.0;

}

OutputTime::~OutputTime(void)
{
    /* It's possible now to do output to the file only in the first process */
     if(rank != 0) {
         /* TODO: do something, when support for Parallel VTK is added */
         return;
     }

     if(this->_base_filename != NULL) {
         delete this->_base_filename;
     }

     if(base_file != NULL) {
         base_file->close();
         delete base_file;
     }

     xprintf(MsgLog, "O.K.\n");
}


void OutputTime::write_all_data(void)
{
    int ierr, rank;
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    ASSERT(ierr == 0, "Error in MPI_Comm_rank.");

    OutputTime *output_time = NULL;

    /* It's possible now to do output to the file only in the first process */
    if(rank != 0) {
        /* TODO: do something, when support for Parallel VTK is added */
        return;
    }

    // Go through all OutputTime objects
    for(std::vector<OutputTime*>::iterator stream_iter = OutputTime::output_streams.begin();
            stream_iter != OutputTime::output_streams.end();
            ++stream_iter)
    {
        // Write data to output stream, when data registered to this output
        // streams were changed
        output_time = (*stream_iter);
        if(output_time->write_time < output_time->time) {
            DBGMSG("Write output to output stream: %s for time: %f\n",
                    (*stream_iter)->name->c_str(),
                    (*stream_iter)->time);
            output_time->write_data();
            // Remember the last time of writing to output stream
            output_time->write_time = output_time->time;
            output_time->current_step++;
        } else {
            DBGMSG("Skipping output stream: %s in time: %f\n",
                    (*stream_iter)->name->c_str(),
                    (*stream_iter)->time);
        }
    }

    /* Free all registered data */
    OutputTime::clear_data();
}


void OutputTime::clear_data(void)
{
	OutputTime *output_time = NULL;

    // Go through all OutputTime objects
    for(std::vector<OutputTime*>::iterator stream_iter = OutputTime::output_streams.begin();
            stream_iter != OutputTime::output_streams.end();
            ++stream_iter)
    {
        output_time = (*stream_iter);
        output_time->node_data.clear();
        output_time->corner_data.clear();
        output_time->elem_data.clear();
    }
}

#define INSTANCE_register_field(spacedim, value) \
	template  void OutputTime::register_data<spacedim, value> \
		(const Input::Record &in_rec, const DiscreteSpace ref_type, Field<spacedim, value> &field);

#define INSTANCE_register_multifield(spacedim, value) \
	template void OutputTime::register_data<spacedim, value> \
		(const Input::Record &in_rec, const DiscreteSpace ref_type, MultiField<spacedim, value> &field);


#define INSTANCE_OutputData(spacedim, value) \
	template class OutputData<value>;


#define INSTANCE_DIM_DEP_VALUES( MACRO, dim_from, dim_to)                      \
		MACRO(dim_from, FieldValue<dim_to>::VectorFixed )                       \
		MACRO(dim_from, FieldValue<dim_to>::TensorFixed )

#define INSTANCE_TO_ALL( MACRO, dim_from) \
		MACRO(dim_from, FieldValue<0>::Enum ) \
		MACRO(dim_from, FieldValue<0>::EnumVector)                \
		MACRO(dim_from, FieldValue<0>::Integer)                \
		MACRO(dim_from, FieldValue<0>::Scalar)                  \
		MACRO(dim_from, FieldValue<0>::Vector)                  \
\
INSTANCE_DIM_DEP_VALUES(MACRO, dim_from, 2) \
INSTANCE_DIM_DEP_VALUES(MACRO, dim_from, 3) \

#define INSTANCE_ALL(MACRO) \
		INSTANCE_TO_ALL( MACRO, 3) \
		INSTANCE_TO_ALL( MACRO, 2)


INSTANCE_ALL( INSTANCE_register_field )
INSTANCE_ALL( INSTANCE_register_multifield )

//INSTANCE_TO_ALL( INSTANCE_OutputData, 0)


//INSTANCE_register_field(3, FieldValue<0>::Scalar)
//INSTANCE_register_multifield(3, FieldValue<0>::Scalar)
//INSTANCE_OutputData(3, FieldValue<0>::Scalar)
