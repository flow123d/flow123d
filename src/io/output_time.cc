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
#include "io/output_data.hh"
#include "io/output_vtk.h"
#include "io/output_msh.h"
#include "mesh/mesh.h"
#include "input/accessors.hh"
#include "system/sys_profiler.hh"


using namespace Input::Type;

Record OutputTime::input_type
    = Record("OutputStream", "Parameters of output.")
    // The stream
    .declare_key("file", FileName::output(), Default::obligatory(),
            "File path to the connected output file.")
            // The format
	.declare_key("format", OutputTime::input_format_type, Default::optional(),
			"Format of output stream and possible parameters.")
	.declare_key("time_step", Double(0.0),
			"Time interval between outputs.\n"
			"Regular grid of output time points starts at the initial time of the equation and ends at the end time which must be specified.\n"
			"The start time and the end time are always added. ")
	.declare_key("time_list", Array(Double(0.0)),
	        Default::read_time("List containing the initial time of the equation. \n You can prescribe an empty list to override this behavior."),
			"Explicit array of output time points (can be combined with 'time_step'.")
	.declare_key("add_input_times", Bool(), Default("false"),
			"Add all input time points of the equation, mentioned in the 'input_fields' list, also as the output points.");


AbstractRecord OutputTime::input_format_type
    = AbstractRecord("OutputTime",
            "Format of output stream and possible parameters.");

OutputDataBase *OutputTime::output_data_by_field_name
		(const std::string &field_name, DiscreteSpace ref_type)
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
        if (data->field_name == field_name)     return data;

    return nullptr;
}

/* Initialize static member of the class */
std::vector<OutputTime*> OutputTime::output_streams;


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
        output_time->format_record_ = *format;
    } else {
        output_time = new OutputVTK(in_rec);
    }

    return output_time;
}

void OutputTime::add_admissible_field_names(const Input::Array &in_array)
{
    vector<Input::FullEnum> field_ids;
    in_array.copy_to(field_ids);

    for (auto it: field_ids) {
        this->output_names.insert(std::pair<std::string, bool>((string)it, true));
    }
}


OutputTime::OutputTime(const Input::Record &in_rec)
: input_record_(in_rec)
{
    ofstream *base_file;

    string fname = in_rec.val<FilePath>("file");

    base_file = new ofstream;

    MPI_Comm_rank(MPI_COMM_WORLD, &this->rank);
    if(this->rank == 0) {
    	base_file->open(fname.c_str());
    	INPUT_CHECK( base_file->is_open() , "Can not open output file: %s\n", fname.c_str() );
    	xprintf(MsgLog, "Writing flow output file: %s ... \n", fname.c_str());
    }

    this->current_step = 0;
    this->_base_file = base_file;
    this->_base_filename = fname;
    this->_data_file = NULL;
    this->_mesh = NULL;
    this->time = -1.0;
    this->write_time = -1.0;
}

OutputTime::~OutputTime(void)
{
    /* It's possible now to do output to the file only in the first process */
     //if(rank != 0) {
     //    /* TODO: do something, when support for Parallel VTK is added */
     //    return;
    // }

     if(this->_base_file != NULL) {
         this->_base_file->close();
         delete this->_base_file;
     }

     xprintf(MsgLog, "O.K.\n");
}


void OutputTime::mark_output_times(const TimeGovernor &tg)
{
	TimeMark::Type output_mark_type = tg.equation_fixed_mark_type() | tg.marks().type_output();

	double time_step;
	if (input_record_.opt_val("time_step", time_step)) {
		tg.add_time_marks_grid(time_step, output_mark_type);
	}

	Input::Array time_list;
	if (input_record_.opt_val("time_list", time_list)) {
		vector<double> list;
		time_list.copy_to(list);
		for( double time : list) tg.marks().add(TimeMark(time, output_mark_type));
	} else {
	    tg.marks().add( TimeMark(tg.init_time(), output_mark_type) );
	}

	bool add_flag;
	if (input_record_.opt_val("add_input_times", add_flag) && add_flag) {
		TimeMark::Type input_mark_type = tg.equation_mark_type() | tg.marks().type_input();
		vector<double> mark_times;
		// can not add marks while iterating through time marks
		for(auto it = tg.marks().begin(input_mark_type); it != tg.marks().end(input_mark_type); ++it)
			mark_times.push_back(it->time());
		for(double time : mark_times)
			tg.marks().add( TimeMark(time, output_mark_type) );

	}

}


void OutputTime::write_time_frame()
{
	START_TIMER("write_time_frame");
    /* TODO: do something, when support for Parallel VTK is added */
    if (this->rank == 0) {
    	// Write data to output stream, when data registered to this output
		// streams were changed
		if(write_time < time) {
			xprintf(MsgLog, "Write output to output stream: %s for time: %f\n",
			        this->_base_filename.c_str(), time);
			write_data();
			// Remember the last time of writing to output stream
			write_time = time;
			current_step++;
		} else {
			xprintf(MsgLog, "Skipping output stream: %s in time: %f\n",
			        this->_base_filename.c_str(), time);
		}
    }
    clear_data();
}

void OutputTime::clear_data(void)
{
	node_data.clear();
    corner_data.clear();
    elem_data.clear();
}



#define INSTANCE_register_field(spacedim, value) \
	template  void OutputTime::register_data<spacedim, value> \
		(const DiscreteSpace ref_type, Field<spacedim, value> &field);

#define INSTANCE_register_multifield(spacedim, value) \
	template void OutputTime::register_data<spacedim, value> \
		(const DiscreteSpace ref_type, MultiField<spacedim, value> &field);


#define INSTANCE_OutputData(spacedim, value) \
	template class OutputData<value>;


#define INSTANCE_DIM_DEP_VALUES( MACRO, dim_from, dim_to) \
		MACRO(dim_from, FieldValue<dim_to>::VectorFixed ) \
		MACRO(dim_from, FieldValue<dim_to>::TensorFixed )

#define INSTANCE_TO_ALL( MACRO, dim_from) \
		MACRO(dim_from, FieldValue<0>::Enum ) \
		MACRO(dim_from, FieldValue<0>::EnumVector) \
		MACRO(dim_from, FieldValue<0>::Integer) \
		MACRO(dim_from, FieldValue<0>::Scalar) \
		MACRO(dim_from, FieldValue<0>::Vector) \
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
