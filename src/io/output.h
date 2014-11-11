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
 * $Id: output.h 2505 2013-09-13 14:52:27Z jiri.hnidek $
 * $Revision: 2505 $
 * $LastChangedBy: jiri.hnidek $
 * $LastChangedDate: 2013-09-13 16:52:27 +0200 (Pá, 13 IX 2013) $
 *
 * @file    output.h
 * @brief   Header: The functions for all outputs.
 *
 *
 * TODO:
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
#include <fstream>
#include <ostream>

#include "system/system.hh"
#include "mesh/mesh.h"

#include "fields/field.hh"
#include "input/accessors.hh"
#include "system/exceptions.hh"
#include "io/output_time.hh"

#include "io/output_data_base.hh"

//class OutputVTK;
//class OutputMSH;



/**
 * \brief This class is used for storing data that are copied from field.
 *
 *
 */
template <class Value>
class OutputData : public OutputDataBase {
public:
	typedef typename Value::element_type ElemType;

	/**
     * \brief Constructor of templated OutputData
     */
	OutputData(const FieldCommon &field,
	        unsigned int size)
	: val_aux(aux)
	{
		this->field_name = field.name();
		this->field_units = field.units();
		this->output_field_name = this->field_name;

		this->n_values=size;
		//val_aux.set_n_comp(10); // just to check that n_elem depends on n_comp

		if (val_aux.NCols_==1) {
			if (val_aux.NRows_==1) {
				this->n_elem_ = scalar;
				this->n_rows = 1;
				this->n_cols = 1;
			} else {
				if (val_aux.NRows_>1) {
					if (val_aux.NRows_ > 3) {
						xprintf(PrgErr, "Do not support output of vectors with fixed size >3. Field: %s\n", this->field_name.c_str());
					} else {
						this->n_elem_ = vector;
						this->n_rows = 3;
						this->n_cols = 1;
					}
				} else {
					THROW(OutputTime::ExcOutputVariableVector() << OutputTime::EI_FieldName(this->field_name));
				}
			}
		} else {
			this->n_elem_ = tensor;
			this->n_rows = 3;
			this->n_cols = 3;
		}


	    data_ = new ElemType[n_values * n_elem_];
	}


    /**
     * \brief Destructor of OutputData
     */
    ~OutputData()
	{
	    delete[] this->data_;
	}


    /**
     * Output data element on given index @p idx. Method for writing data to output stream
     *
     * TODO: should at least output whole output value at once, since storage format should be hidden.
     * TODO: should output whole array at once, otherwise this could be performance bottleneck.
     * TODO: indicate if the tensor data are output in column-first or raw-first order
     *       and possibly implement transposition. Set such property for individual file formats.
     *       Class OutputData stores always in raw-first order.
     */
    void print(ostream &out_stream, unsigned int idx) override
    {
        ASSERT_LESS(idx, this->n_values);
        ElemType *ptr_begin = data_ + n_elem_ * idx;
        for(ElemType *ptr = ptr_begin; ptr < ptr_begin + n_elem_; ptr++ )
        	out_stream << *ptr << " ";
    }

    /**
     * Store data element of given data value under given index.
     */
    void store_value(unsigned int idx, const Value& value) {
    	operate(idx, value,  [](ElemType& raw, ElemType val) {raw=val;});
    };

    /**
     * Add value to given index
     */
    void add(unsigned int idx, const Value& value) {
    	operate(idx, value,   [](ElemType& raw, ElemType val) {raw+=val;});
    };

    void zero(unsigned int idx) {
    	operate(idx, val_aux, 	[](ElemType& raw, ElemType val) {raw=0;});
    };

    void normalize(unsigned int idx, unsigned int divisor) {
    	operate(idx, val_aux, 	[divisor](ElemType& raw, ElemType val) {raw/=divisor;});
    };



private:


    template <class Func>
    void operate(unsigned int idx, const Value &val, const Func& func) {
    	ASSERT_LESS(idx, this->n_values);
    	ElemType *ptr = data_ + idx*n_elem_;
        for(unsigned int i_row=0; i_row < this->n_rows; i_row++)
        	for(unsigned int i_col=0; i_col < this->n_cols; i_col++)
        	{
        		if (i_row < val.n_rows() && i_col < val.n_cols())
        			func(*ptr, val(i_row, i_col));
        		else
        			func(*ptr, 0);
        		ptr++;
        	}
    };





    /**
     * Computed data values for output stored as continuous buffer of their data elements.
     * One data value has @p n_elem data elements (of type double, int or unsigned int).
     */
    ElemType *data_;

    /// auxiliary value
    typename Value::return_type aux;
    // auxiliary field value envelope over @p aux
    Value val_aux;

    // Number of rows and cols in stored data element, valid values are
    // (1,1) for scalar; (3,1) for vectors; (3,3) for tensors
    unsigned int n_rows, n_cols;



};



/**************************************************************************************************************
 * OutputTime implementation
 */

template<int spacedim, class Value>
void OutputTime::register_data(const DiscreteSpace type,
        MultiField<spacedim, Value> &multi_field)
{
	if (output_names.find(multi_field.name()) != output_names.end()) {
		if (output_names[multi_field.name()] == true)
			for (unsigned long index=0; index < multi_field.size(); index++)
				this->compute_field_data(type, multi_field[index] );

		return;
	}
	else
	{
		// We are trying to output field that is not recognized by the stream.
		//DBGMSG("Internal error: Output stream %s does not support field %s.\n", name.c_str(), multi_field.name().c_str());
	}
}


template<int spacedim, class Value>
void OutputTime::register_data(const DiscreteSpace ref_type,
        Field<spacedim, Value> &field_ref)
{
	if (output_names.find(field_ref.name()) != output_names.end()) {
		if (output_names[field_ref.name()] == true)
			this->compute_field_data(ref_type, field_ref);

		return;
	}
	else
	{
		// We are trying to output field that is not recognized by the stream.
		//DBGMSG("Internal error: Output stream %s does not support field %s.\n", name.c_str(), field_ref.name().c_str());
	}
}


template<int spacedim, class Value>
void OutputTime::compute_field_data(DiscreteSpace space_type, Field<spacedim, Value> &field)
{

    /* It's possible now to do output to the file only in the first process */
    if( this->rank != 0) {
        /* TODO: do something, when support for Parallel VTK is added */
        return;
    }


    // TODO: remove const_cast after resolving problems with const Mesh.
    mesh = const_cast<Mesh *>(field.mesh());
    ASSERT(mesh, "Null mesh pointer.\n");

    // get possibly existing data for the same field, check both name and type
    OutputDataBase *data = output_data_by_field_name(field.name(), space_type);
    OutputData<Value> *output_data = dynamic_cast<OutputData<Value> *>(data);

    if (!output_data) {
        switch(space_type) {
        case NODE_DATA:
        	output_data = new OutputData<Value>(field, mesh->n_nodes());
            node_data.push_back(output_data);
            break;
        case CORNER_DATA: {
            unsigned int n_corners = 0;
            FOR_ELEMENTS(mesh, ele)
                n_corners += ele->n_nodes();
        	output_data = new OutputData<Value>(field, n_corners );
            corner_data.push_back(output_data);
        }
        break;
        case ELEM_DATA:
        	output_data = new OutputData<Value>(field, mesh->n_elements() );
            elem_data.push_back(output_data);
            break;
        }
    }

    unsigned int i_node;

    /* Copy data to array */
    switch(space_type) {
    case NODE_DATA: {
    	// set output data to zero
    	vector<unsigned int> count(output_data->n_values, 0);
    	for(unsigned int idx=0; idx < output_data->n_values; idx++)
    		output_data->zero(idx);

    	// sum values
        FOR_ELEMENTS(mesh, ele) {
            FOR_ELEMENT_NODES(ele, i_node) {
                Node * node = ele->node[i_node];
                unsigned int ele_index = ele.index();
                unsigned int node_index = mesh->node_vector.index(ele->node[i_node]);

				const Value &node_value =
						Value( const_cast<typename Value::return_type &>(
								field.value(node->point(), ElementAccessor<spacedim>(mesh, ele_index,false)) ));
				output_data->add(node_index, node_value);
				count[node_index]++;

            }
        }

        // Compute mean values at nodes
    	for(unsigned int idx=0; idx < output_data->n_values; idx++)
    		output_data->normalize(idx, count[idx]);
    }
    break;
    case CORNER_DATA: {
    	unsigned int corner_index=0;
        FOR_ELEMENTS(mesh, ele) {
            FOR_ELEMENT_NODES(ele, i_node) {
                Node * node = ele->node[i_node];
                unsigned int ele_index = ele.index();

				const Value &node_value =
						Value( const_cast<typename Value::return_type &>(
								field.value(node->point(), ElementAccessor<spacedim>(mesh, ele_index,false)) ));
                output_data->store_value(corner_index,  node_value);
                corner_index++;
            }
        }
    }
    break;
    case ELEM_DATA: {
        FOR_ELEMENTS(mesh, ele) {
            unsigned int ele_index = ele.index();
			const Value &ele_value =
					Value( const_cast<typename Value::return_type &>(
							field.value(ele->centre(), ElementAccessor<spacedim>(mesh, ele_index,false)) ));
			//std::cout << ele_index << " ele:" << typename Value::return_type(ele_value) << std::endl;
            output_data->store_value(ele_index,  ele_value);
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
