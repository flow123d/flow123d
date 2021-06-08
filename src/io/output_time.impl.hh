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

#include "system/system.hh"

#include "system/exceptions.hh"
#include "io/output_time.hh"
#include "io/output_mesh.hh"
#include "io/element_data_cache.hh"
#include "io/output_element.hh"




/**************************************************************************************************************
 * OutputTime implementation
 */

template <typename T>
OutputTime::OutputDataPtr OutputTime::prepare_compute_data(std::string field_name,
    DiscreteSpace space_type, unsigned int n_rows, unsigned int n_cols)
{
    // get possibly existing data for the same field, check both name and type
    unsigned int size;
    switch (space_type) {
        case NODE_DATA:
        case CORNER_DATA:
        {
            auto &offset_vec = *( output_mesh_->offsets_->get_component_data(0).get() );
            size = offset_vec[offset_vec.size()-1];
            break;
        }
        case ELEM_DATA:
        case NATIVE_DATA:
            size = output_mesh_->offsets_->n_values()-1;
            break;
        default:
            ASSERT(false).error("Should not happen.");
            break;
    }

    /**
     * supposing that everything is output at the first step !
     * - if(current_step == 0) push_back all the non-existing fields shared<ElementDataCache<T>>
     * - else just set the correct shared<ElementDataCache<T>>
     * - when clearing output_data_vec, set all to shared<DummyElementCache>
     * - this can be done in many equations sharing the same OutputTime stream
     * 
    */
    auto &od_vec=this->output_data_vec_[space_type];
    auto it=std::find_if(od_vec.begin(), od_vec.end(),
            [&field_name](OutputDataPtr ptr) { return (ptr->field_input_name() == field_name); });
    
    if(current_step == 0 && it == od_vec.end() ) {
        // DebugOut() << "OutputTime::prepare_compute_data: PUSH BACK " << field_name << "\n";
        od_vec.push_back( std::make_shared< ElementDataCache<T> >(field_name, n_rows*n_cols, size) );
        it=--od_vec.end();
    }
    else{
        ASSERT(it != od_vec.end()).error("Try to add non-existing output field after first step.");
        *it = std::make_shared< ElementDataCache<T> >(field_name, n_rows*n_cols, size);
    }

    return *it;
}

#endif
