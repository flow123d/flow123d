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
 * @file    output.cc
 * @brief   The functions for all outputs (methods of classes: Output and OutputTime).
 */

#include <string>

#include "system/sys_profiler.hh"
#include "mesh/mesh.h"
#include "input/accessors.hh"
#include "output_time.impl.hh"
#include "output_vtk.hh"
#include "output_msh.hh"
#include "output_mesh.hh"
#include "io/output_time_set.hh"
#include "io/observe.hh"


FLOW123D_FORCE_LINK_IN_PARENT(vtk)
FLOW123D_FORCE_LINK_IN_PARENT(gmsh)


namespace IT = Input::Type;

const IT::Record & OutputTime::get_input_type() {
    return IT::Record("OutputStream", "Configuration of the spatial output of a single balance equation.")
		// The stream
		.declare_key("file", IT::FileName::output(), IT::Default::read_time("Name of the equation associated with the output stream."),
				"File path to the connected output file.")
				// The format
		.declare_key("format", OutputTime::get_input_format_type(), IT::Default("{}"),
				"Format of output stream and possible parameters.")
		.declare_key("times", OutputTimeSet::get_input_type(), IT::Default::optional(),
		        "Output times used for equations without is own output times key.")
        .declare_key("output_mesh", OutputMeshBase::get_input_type(), IT::Default::optional(),
                "Output mesh record enables output on a refined mesh [EXPERIMENTAL, VTK only]."
                "Sofar refinement is performed only in discontinous sense."
                "Therefore only corner and element data can be written on refined output mesh.
                "Node data are to be transformed to corner data, native data cannot be written."
                "Do not include any node or native data in output fields.")
        .declare_key("precision", IT::Integer(0), IT::Default("5"),
                "The number of decimal digits used in output of floating point values.")
        .declare_key("observe_points", IT::Array(ObservePoint::get_input_type()), IT::Default("[]"),
                "Array of observe points.")
		.close();
}


IT::Abstract & OutputTime::get_input_format_type() {
	return IT::Abstract("OutputTime", "Format of output stream and possible parameters.")
	    .allow_auto_conversion("vtk")
		.close();
}


OutputTime::OutputTime()
: current_step(0),
  time(-1.0),
  write_time(-1.0),
  _mesh(nullptr)
{
    MPI_Comm_rank(MPI_COMM_WORLD, &this->rank);
}



void OutputTime::init_from_input(const std::string &equation_name, Mesh &mesh, const Input::Record &in_rec)
{
    _mesh = &mesh;

    input_record_ = in_rec;
    equation_name_ = equation_name;

    // Read output base file name
    // TODO: remove dummy ".xyz" extension after merge with DF
    FilePath output_file_path(equation_name+"_fields", FilePath::output_file);
    input_record_.opt_val("file", output_file_path);
    this->_base_filename = output_file_path;
}



OutputTime::~OutputTime(void)
{
    /* It's possible now to do output to the file only in the first process */
     //if(rank != 0) {
     //    /* TODO: do something, when support for Parallel VTK is added */
     //    return;
    // }

    if (this->_base_file.is_open()) this->_base_file.close();

    LogOut() << "O.K.";
}


Input::Iterator<Input::Array> OutputTime::get_time_set_array() {
    return input_record_.find<Input::Array>("times");
}


Input::Iterator<Input::Record> OutputTime::get_output_mesh_record() {
    return input_record_.find<Input::Record>("output_mesh");
}


std::shared_ptr<OutputMeshBase> OutputTime::create_output_mesh_ptr(bool init_input, bool discont) {
	if (discont) {
		if (init_input) output_mesh_discont_ = std::make_shared<OutputMeshDiscontinuous>(*_mesh, *this->get_output_mesh_record());
		else output_mesh_discont_ = std::make_shared<OutputMeshDiscontinuous>(*_mesh);
		return output_mesh_discont_;
	} else {
		if (init_input) output_mesh_ = std::make_shared<OutputMesh>(*_mesh, *this->get_output_mesh_record());
		else output_mesh_ = std::make_shared<OutputMesh>(*_mesh);
		return output_mesh_;
	}
}


std::shared_ptr<OutputMeshBase> OutputTime::get_output_mesh_ptr(bool discont) {
	if (discont) {
        compute_discontinuous_output_mesh();
		return output_mesh_discont_;
	} else {
		return output_mesh_;
	}
}


void OutputTime::update_time(double field_time) {
	if (this->time < field_time) {
		this->time = field_time;
	}
}


void OutputTime::compute_discontinuous_output_mesh()
{
    ASSERT_PTR(output_mesh_).error("Create output mesh first!");
    if( ! output_mesh_discont_->is_refined())
        output_mesh_discont_->create_mesh(output_mesh_);
}



void OutputTime::fix_main_file_extension(std::string extension)
{
    if(extension.compare( this->_base_filename.extension() ) != 0) {
        string old_name = (string)this->_base_filename;
        std::vector<string> path = {this->_base_filename.parent_path(), this->_base_filename.stem() + extension};
        this->_base_filename = FilePath(
                path,
                FilePath::output_file);
        WarningOut() << "Renaming output file: " << old_name << " to " << this->_base_filename;

    }
}


/* Initialize static member of the class */
//std::vector<OutputTime*> OutputTime::output_streams;


/*
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
    */


std::shared_ptr<OutputTime> OutputTime::create_output_stream(const std::string &equation_name, Mesh &mesh, const Input::Record &in_rec)
{

    Input::AbstractRecord format = Input::Record(in_rec).val<Input::AbstractRecord>("format");
    std::shared_ptr<OutputTime> output_time = format.factory< OutputTime >();
    output_time->init_from_input(equation_name, mesh, in_rec);

    return output_time;
}





void OutputTime::write_time_frame()
{
	START_TIMER("OutputTime::write_time_frame");
    /* TODO: do something, when support for Parallel VTK is added */
    if (observe_)
        observe_->output_time_frame(time);

    if (this->rank == 0) {

    	// Write data to output stream, when data registered to this output
		// streams were changed
		if(write_time < time) {

			LogOut() << "Write output to output stream: " << this->_base_filename << " for time: " << time;
			write_data();
			// Remember the last time of writing to output stream
			write_time = time;
			current_step++;
            
            // invalidate output meshes after the time frame written
            output_mesh_.reset();
            output_mesh_discont_.reset();
		} else {
			LogOut() << "Skipping output stream: " << this->_base_filename << " in time: " << time;
		}
    }
    clear_data();
}

std::shared_ptr<Observe> OutputTime::observe()
{
    ASSERT_PTR(_mesh);
    // create observe object at first call
    if (! observe_) {
        auto observe_points = input_record_.val<Input::Array>("observe_points");
        unsigned int precision = input_record_.val<unsigned int>("precision");
        observe_ = std::make_shared<Observe>(this->equation_name_, *_mesh, observe_points, precision);
    }
    return observe_;
}


void OutputTime::clear_data(void)
{
    for(auto &map : output_data_vec_)  map.clear();
}



// explicit instantiation of template methods
#define OUTPUT_PREPARE_COMPUTE_DATA(TYPE) \
template ElementDataCache<TYPE> & OutputTime::prepare_compute_data<TYPE>(std::string field_name, DiscreteSpace space_type, \
		unsigned int n_rows, unsigned int n_cols)

OUTPUT_PREPARE_COMPUTE_DATA(int);
OUTPUT_PREPARE_COMPUTE_DATA(unsigned int);
OUTPUT_PREPARE_COMPUTE_DATA(double);

