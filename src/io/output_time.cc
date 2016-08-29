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
#include <fields/field_set.hh>


FLOW123D_FORCE_LINK_IN_PARENT(vtk)
FLOW123D_FORCE_LINK_IN_PARENT(gmsh)


namespace IT = Input::Type;

const IT::Record & OutputTime::get_input_type() {
    return IT::Record("OutputStream", "Parameters of output.")
		// The stream
		.declare_key("file", IT::FileName::output(), IT::Default::read_time("Name of the equation associated with the output stream."),
				"File path to the connected output file.")
				// The format
		.declare_key("format", OutputTime::get_input_format_type(), IT::Default("{}"),
				"Format of output stream and possible parameters.")
		.declare_key("times", OutputTimeSet::get_input_type(), IT::Default::optional(),
		        "Output times used for equations without is own output times key.")
        .declare_key("output_mesh", OutputMeshBase::get_input_type(), IT::Default::optional(),
                "Output mesh record enables output on a refined mesh.")
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


void OutputTime::make_output_mesh(FieldSet &output_fields)
{

    // make observe points if not already done
    observe();

    // already computed
    if(output_mesh_) return;

    // Read optional error control field name
    auto it = input_record_.find<Input::Record>("output_mesh");
    
    if(enable_refinement_) {
        if(it) {
            output_mesh_ = std::make_shared<OutputMesh>(*_mesh, *it);
            output_mesh_discont_ = std::make_shared<OutputMeshDiscontinuous>(*_mesh, *it);
            output_mesh_->select_error_control_field(output_fields);
            output_mesh_discont_->select_error_control_field(output_fields);
            
            output_mesh_->create_refined_mesh();
            return;
        }
    }
    else
    {
        // skip creation of output mesh (use computational one)
        if(it)
        	WarningOut() << "Ignoring output mesh record.\n Output in GMSH format available only on computational mesh!";
    }
    
    
    output_mesh_ = std::make_shared<OutputMesh>(*_mesh);
    output_mesh_discont_ = std::make_shared<OutputMeshDiscontinuous>(*_mesh);
    
    output_mesh_->create_identical_mesh();
}


void OutputTime::compute_discontinuous_output_mesh()
{
    ASSERT_PTR(output_mesh_).error("Create output mesh first!");
    output_mesh_discont_->create_mesh(output_mesh_);
}



void OutputTime::fix_main_file_extension(std::string extension)
{
    if(extension.compare( this->_base_filename.extension() ) != 0) {
        string new_name = (string)this->_base_filename + extension;
        WarningOut() << "Renaming output file: " << this->_base_filename << " to " << new_name;
        this->_base_filename = new_name;
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
		MACRO(dim_from, FieldValue<0>::Integer) \
		MACRO(dim_from, FieldValue<0>::Scalar) \
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
