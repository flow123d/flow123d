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
 * @file    output.h
 * @brief   Header: The functions for all outputs.
 * @todo
 * - remove Output, keep OutputTime only (done)
 * - remove parameter mesh from static method OutputTime::output_stream (done)
 * - move initialization of streams from hc_expolicit_sequantial to
 *     Aplication::Aplication() constructor (done)
 * - OutputTime::register_XXX_data - should accept iterator to output record of particular equation, ask for presence of the key
 *   that has same name as the name of the quantity to output, extract the string with stream name from this key, find the stream
 *   and perform output.
 *
 *   on input:
 *
 *   { // darcy flow
 *      output = {
 *          pressure_nodes="nodal_data",
 *          pressure_elements="el_data"
 *      }
 *   }
 *
 *   output_streams=[
 *      {name="nodal_data", ... },
 *      {name="el_data", ... }
 *   ]
 *
 *   in code:
 *
 *   Input::Record out_rec = in_rec.val<Input::Record>("output");
 *   OutputTime::register_node_data(mesh_, "pressure_nodes", "L", out_rec, node_pressure);
 *   OutputTime::register_elem_data(mesh_, "pressure_elements", "L", out_rec, ele_pressure);
 *   ...
 *
 * - use exceptions instead of returning result, see declaration of exceptions through DECLARE_EXCEPTION macro
 * - move write_data from equations into coupling, write all streams
 *
 * =======================
 * - Is it still necessary to split output into registration and write the data?
 *   Could we perform it at once? ... No, it doesn't make any sense.
 * - Support for output of corner data into GMSH format (ElementNodeData section)
 *
 */

#ifndef OUTPUT_H
#define OUTPUT_H

#include <vector>
#include <string>
#include <ostream>

#include "system/system.hh"
#include "mesh/mesh.h"

#include "fields/field.hh"
#include "fields/multi_field.hh"
#include "system/exceptions.hh"
#include "io/output_time.hh"

#include "io/output_data_base.hh"
#include "output_mesh.hh"
#include "io/output_data.hh"
#include "output_element.hh"





/**************************************************************************************************************
 * OutputTime implementation
 */

template<int spacedim, class Value>
void OutputTime::register_data(const DiscreteSpace type,
        MultiField<spacedim, Value> &multi_field)
{
	OLD_ASSERT_LESS(type, N_DISCRETE_SPACES);

    DiscreteSpaceFlags flags = 1 << type;
    for (unsigned long index=0; index < multi_field.size(); index++)
        for(unsigned int ids=0; ids < N_DISCRETE_SPACES; ids++)
            if (flags & (1 << ids))
                    this->compute_field_data( DiscreteSpace(ids), multi_field[index] );

}


template<int spacedim, class Value>
void OutputTime::register_data(const DiscreteSpace type,
        Field<spacedim, Value> &field_ref)
{
	OLD_ASSERT_LESS(type, N_DISCRETE_SPACES);
    
	DiscreteSpaceFlags flags = 1 << type;
    for(unsigned int ids=0; ids < N_DISCRETE_SPACES; ids++)
        if (flags & (1 << ids))
            this->compute_field_data( DiscreteSpace(ids), field_ref);
}


template<int spacedim, class Value>
void OutputTime::compute_field_data(DiscreteSpace space_type, Field<spacedim, Value> &field)
{
    /* It's possible now to do output to the file only in the first process */
    if( this->rank != 0) {
        /* TODO: do something, when support for Parallel VTK is added */
        return;
    }

    if(space_type == CORNER_DATA)
        compute_discontinuous_output_mesh();
    
    // get possibly existing data for the same field, check both name and type
    std::vector<unsigned int> size(N_DISCRETE_SPACES);
    size[NODE_DATA] = output_mesh_->n_nodes();
    size[ELEM_DATA] = output_mesh_->n_elements();
    size[CORNER_DATA] = output_mesh_discont_->n_nodes();

    auto &od_vec=this->output_data_vec_[space_type];
    auto it=std::find_if(od_vec.begin(), od_vec.end(),
            [&field](OutputDataPtr ptr) { return (ptr->field_name ==  field.name()); });
    if ( it == od_vec.end() ) {
        od_vec.push_back( std::make_shared< OutputData<Value> >(field, size[space_type]) );
        it=--od_vec.end();
    }
    OutputData<Value> &output_data = dynamic_cast<OutputData<Value> &>(*(*it));


    /* Copy data to array */
    switch(space_type) {
    case NODE_DATA: {
        //DBGMSG("compute field NODE data\n");
        // set output data to zero
        vector<unsigned int> count(output_data.n_values, 0);
        for(unsigned int idx=0; idx < output_data.n_values; idx++)
            output_data.zero(idx);

        for(const auto & ele : *output_mesh_)
        {
            std::vector<Space<3>::Point> vertices = ele.vertex_list();
            for(unsigned int i=0; i < ele.n_nodes(); i++)
            {
                unsigned int node_index = ele.node_index(i);
                const Value &node_value =
                        Value( const_cast<typename Value::return_type &>(
                                field.value(vertices[i],
                                            ElementAccessor<spacedim>(ele.orig_mesh(), ele.orig_element_idx(),false) ))
                             );
                output_data.add(node_index, node_value);
                count[node_index]++;
            }
        }
        
        // Compute mean values at nodes
        for(unsigned int idx=0; idx < output_data.n_values; idx++)
            output_data.normalize(idx, count[idx]);
    }
    break;
    case CORNER_DATA: {
        for(const auto & ele : *output_mesh_discont_)
        {
            DBGMSG("ele %d\n",ele.idx());
            std::vector<Space<3>::Point> vertices = ele.vertex_list();
            for(unsigned int i=0; i < ele.n_nodes(); i++)
            {
                unsigned int node_index = ele.node_index(i);
                const Value &node_value =
                        Value( const_cast<typename Value::return_type &>(
                                field.value(vertices[i],
                                            ElementAccessor<spacedim>(ele.orig_mesh(), ele.orig_element_idx(),false) ))
                             );
                output_data.store_value(node_index,  node_value);
            }
        }
        DBGMSG("compute field CORNER data END\n");
    }
    break;
    case ELEM_DATA: {
        for(const auto & ele : *output_mesh_)
        {
            unsigned int ele_index = ele.idx();
            const Value &ele_value =
                        Value( const_cast<typename Value::return_type &>(
                                field.value(ele.centre(),
                                            ElementAccessor<spacedim>(ele.orig_mesh(), ele.orig_element_idx(),false))
                                                                        )
                             );
                output_data.store_value(ele_index,  ele_value);
        }
    }
    break;
    }

    /* Set the last time */
    if(this->time < field.time()) {
        this->time = field.time();
    }
}


#endif
