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
    stringstream default_prec;
    default_prec << std::numeric_limits<double>::max_digits10;
    return IT::Record("OutputStream", "Configuration of the spatial output of a single balance equation.")
		// The stream
		.declare_key("file", IT::FileName::output(), IT::Default::read_time("Name of the equation associated with the output stream."),
				"File path to the connected output file.")
				// The format
		.declare_key("format", OutputTime::get_input_format_type(), IT::Default("{}"),
				"File format of the output stream and possible parameters.")
		.declare_key("times", OutputTimeSet::get_input_type(), IT::Default::optional(),
		        "Output times used for fields that do not have their own output times defined.")
        .declare_key("output_mesh", OutputMeshBase::get_input_type(), IT::Default::optional(),
                "Output mesh record enables output on a refined mesh [EXPERIMENTAL, VTK only]."
                "Sofar refinement is performed only in discontinous sense."
                "Therefore only corner and element data can be written on refined output mesh."
                "Node data are to be transformed to corner data, native data cannot be written."
                "Do not include any node or native data in output fields.")
        .declare_key("precision", IT::Integer(0), IT::Default(default_prec.str()),
                "The number of decimal digits used in output of floating point values.\n"
                "Default is 17 decimal digits which are necessary to reproduce double values exactly after write-read cycle.")
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
  parallel_(false)
{
    MPI_Comm_rank(MPI_COMM_WORLD, &this->rank_);
    MPI_Comm_size(MPI_COMM_WORLD, &this->n_proc_);
}



void OutputTime::init_from_input(const std::string &equation_name, const Input::Record &in_rec, std::string unit_str)
{

    input_record_ = in_rec;
    equation_name_ = equation_name;
    unit_string_ = unit_str;

    // Read output base file name
    // TODO: remove dummy ".xyz" extension after merge with DF
    FilePath output_file_path(equation_name+"_fields", FilePath::output_file);
    input_record_.opt_val("file", output_file_path);
    this->precision_ = input_record_.val<int>("precision");
    this->_base_filename = output_file_path;
}


void OutputTime::set_stream_precision(std::ofstream &stream)
{
    //stream.setf(std::ios::scientific);
    stream.precision(this->precision_);
}


OutputTime::~OutputTime(void)
{
    /* It's possible now to do output to the file only in the first process */
     //if(rank_ != 0) {
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


void OutputTime::set_output_data_caches(std::shared_ptr<OutputMeshBase> mesh_ptr) {
    this->nodes_ = mesh_ptr->get_master_mesh()->nodes_;
    this->connectivity_ = mesh_ptr->get_master_mesh()->connectivity_;
    this->offsets_ = mesh_ptr->get_master_mesh()->offsets_;
    output_mesh_ = mesh_ptr;
}


std::shared_ptr<OutputMeshBase> OutputTime::get_output_mesh_ptr() {
	return output_mesh_;
}


void OutputTime::update_time(double field_time) {
	if (this->time < field_time) {
		this->time = field_time;
	}
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


std::shared_ptr<OutputTime> OutputTime::create_output_stream(const std::string &equation_name, const Input::Record &in_rec, std::string unit_str)
{

    Input::AbstractRecord format = Input::Record(in_rec).val<Input::AbstractRecord>("format");
    std::shared_ptr<OutputTime> output_time = format.factory< OutputTime >();
    output_time->init_from_input(equation_name, in_rec, unit_str);

    return output_time;
}





void OutputTime::write_time_frame()
{
	START_TIMER("OutputTime::write_time_frame");
    if (observe_)
        observe_->output_time_frame( write_time < time );

    // Write data to output stream, when data registered to this output
    // streams were changed
    if(write_time < time) {

    	if (this->rank_ == 0 || this->parallel_) // for serial output write log only one (same output file on all processes)
    	    LogOut() << "Write output to output stream: " << this->_base_filename << " for time: " << time;
    	gather_output_data();
        write_data();
        // Remember the last time of writing to output stream
        write_time = time;
        current_step++;
            
        // invalidate output data caches after the time frame written
        // TODO we need invalidate pointers only in special cases (e. g. refining of mesh)
        /*output_mesh_.reset();
        this->nodes_.reset();
        this->connectivity_.reset();
        this->offsets_.reset();*/
    } else {
    	if (this->rank_ == 0 || this->parallel_) // for serial output write log only one (same output file on all processes)
    	    LogOut() << "Skipping output stream: " << this->_base_filename << " in time: " << time;
    }
    clear_data();
}

std::shared_ptr<Observe> OutputTime::observe(Mesh *mesh)
{
    // create observe object at first call
    if (! observe_) {
        auto observe_points = input_record_.val<Input::Array>("observe_points");
        unsigned int precision = input_record_.val<unsigned int>("precision");
        observe_ = std::make_shared<Observe>(this->equation_name_, *mesh, observe_points, precision, this->unit_string_);
    }
    return observe_;
}


void OutputTime::clear_data(void)
{
    for(auto &map : output_data_vec_)  map.clear();
}


int OutputTime::get_parallel_current_step()
{
	if (parallel_) return n_proc_*current_step+rank_;
	else return current_step;
}


void OutputTime::add_dummy_fields()
{}


void OutputTime::gather_output_data(void)
{
    /* for serial output call gather of all data sets */
    if ( !parallel_ ) {
    	auto &offset_vec = *( output_mesh_->offsets_->get_component_data(0).get() );

    	auto &node_data_map = this->output_data_vec_[NODE_DATA];
        for(unsigned int i=0; i<node_data_map.size(); ++i) {
            auto elem_node_cache = node_data_map[i]->element_node_cache_fixed_size(offset_vec);
            auto serial_fix_data_cache = elem_node_cache->gather(output_mesh_->el_ds_, output_mesh_->el_4_loc_);
            if (rank_==0) {
            	auto &master_offset_vec = *( this->offsets_->get_component_data(0).get() );
            	auto serial_data_cache = serial_fix_data_cache->element_node_cache_optimize_size(master_offset_vec);
            	auto &master_conn_vec = *( this->connectivity_->get_component_data(0).get() );
            	node_data_map[i] = serial_data_cache->compute_node_data(master_conn_vec, this->nodes_->n_values());
            }
        }

        auto &corner_data_map = this->output_data_vec_[CORNER_DATA];
        for(unsigned int i=0; i<corner_data_map.size(); ++i) {
            auto elem_node_cache = corner_data_map[i]->element_node_cache_fixed_size(offset_vec);
            auto serial_fix_data_cache = elem_node_cache->gather(output_mesh_->el_ds_, output_mesh_->el_4_loc_);
            if (rank_==0) {
                auto &master_offset_vec = *( this->offsets_->get_component_data(0).get() );
                corner_data_map[i] = serial_fix_data_cache->element_node_cache_optimize_size(master_offset_vec);
            }
        }

    	auto &elm_data_map = this->output_data_vec_[ELEM_DATA];
    	for(unsigned int i=0; i<elm_data_map.size(); ++i) {
    	    auto serial_data = elm_data_map[i]->gather(output_mesh_->el_ds_, output_mesh_->el_4_loc_);
    	    if (rank_==0) elm_data_map[i] = serial_data;
    	}
    	auto &native_data_map = this->output_data_vec_[NATIVE_DATA];
    	for(unsigned int i=0; i<native_data_map.size(); ++i) {
    	    auto serial_data = native_data_map[i]->gather(output_mesh_->el_ds_, output_mesh_->el_4_loc_);
    	    if (rank_==0) {
    	    	auto hash = native_data_map[i]->dof_handler_hash();
    	        native_data_map[i] = serial_data;
    	        (native_data_map[i])->set_dof_handler_hash(hash);
    	    }
    	}
    } else {
        /* Parallel output needs compute node data (average values) */
    	auto &conn_vec = *( output_mesh_->connectivity_->get_component_data(0).get() );
    	auto &node_data_map = this->output_data_vec_[NODE_DATA];
    	for(unsigned int i=0; i<node_data_map.size(); ++i) {
    		node_data_map[i] = node_data_map[i]->compute_node_data(conn_vec, this->nodes_->n_values());
    	}
    }
}



// explicit instantiation of template methods
#define OUTPUT_PREPARE_COMPUTE_DATA(TYPE) \
template ElementDataCache<TYPE> & OutputTime::prepare_compute_data<TYPE>(std::string field_name, DiscreteSpace space_type, \
		unsigned int n_rows, unsigned int n_cols)

OUTPUT_PREPARE_COMPUTE_DATA(int);
OUTPUT_PREPARE_COMPUTE_DATA(unsigned int);
OUTPUT_PREPARE_COMPUTE_DATA(double);

