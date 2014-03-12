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
#include <typeinfo>
#include <mpi.h>
#include <boost/any.hpp>
#include <assert.h>

#include "system/xio.h"
#include "mesh/mesh.h"

#include "fields/field.hh"
#include "input/accessors.hh"
#include "system/exceptions.hh"

class OutputVTK;
class OutputMSH;

/**
 * \brief This method is generic parent class for templated OutputData
 */
class OutputDataBase {
public:
    OutputDataBase() {
        this->field = NULL;
        this->items_count = 0;
        this->vector_items_count = 0;
    };
    virtual ~OutputDataBase() {};

    /**
     *
     */
    FieldCommonBase *field;
    /**
     *
     */
    int items_count;

    /**
     *
     */
    int vector_items_count;

    /**
     *
     */
    virtual void print(ostream &out_stream, unsigned int idx) = 0;
protected:

};

/**
 * \brief This class is used for storing data that are copied from field.
 */
template <class T>
class OutputData : public OutputDataBase {
public:
    /**
     * \brief Constructor of templated OutputData
     */
    OutputData(FieldCommonBase *field, int items_count, int vector_items_count);

    /**
     * \brief Destructor of OutputData
     */
    ~OutputData();

    /**
     * Method for writing data to output stream
     */
    void print(ostream &out_stream, unsigned int idx) {
        assert(idx < this->items_count);
        out_stream << this->_data[idx];
    }

    /**
     * Overloaded operator []
     */
    T& operator[] (std::size_t idx) {
        assert(idx < this->items_count);
        return this->_data[idx];
    };

private:
    /**
     * Array of templated data
     */
    T *_data;
};

template <class T>
OutputData<T>::OutputData(FieldCommonBase *field,
        int items_count,
        int vector_items_count)
{
    this->field = field;
    this->items_count = items_count;
    this->_data = new T[vector_items_count * items_count];
    this->vector_items_count = vector_items_count;
}

template <class T>
OutputData<T>::~OutputData()
{
    delete[] this->_data;
}

/**
 * \brief The class for outputing data during time.
 *
 * This class is descendant of Output class. This class is used for outputing
 * data varying in time. Own output to specific file formats is done at other
 * places to. See output_vtk.cc and output_msh.cc.
 */
class OutputTime {

public:
    /**
     * Declaration of exceptions
     */
    //TYPEDEF_ERR_INFO(EI_StreamName, const string);
    //DECLARE_INPUT_EXCEPTION( ExcStreamName, << "Not valid stream name " << EI_StreamName::qval << "\n");

    /**
     * \brief Constructor of OutputTime object. It opens base file for writing.
     *
     * \param[in] in_rec The reference on the input record
     */
    OutputTime(const Input::Record &in_rec);

    /**
     * \brief Destructor of OutputTime. It doesn't do anything, because all
     * necessary destructors will be called in destructor of Output
     */
    virtual ~OutputTime();

    vector<OutputDataBase*>    node_data;
    vector<OutputDataBase*>    corner_data;
    vector<OutputDataBase*>    elem_data;

    double          time;               ///< The newest time of registered data

    double          write_time;         ///< The last time, when data was wrote to this stream

    ofstream& get_base_file(void) { return *base_file; };

    ofstream& get_data_file(void) { return *data_file; };

    string& get_data_filename(void) { return *data_filename; };

    Mesh *get_mesh(void) { return mesh; };

    unsigned int get_corner_count(void) {
        unsigned int li, count = 0;
        FOR_ELEMENTS(this->mesh, ele) {
            FOR_ELEMENT_NODES(ele, li) {
                count++;
            }
        }
        return count;
    }

    void set_data_file(ofstream *_data_file) { data_file = _data_file; };

    /**
     * Enumeration of file formats supported by Flow123d
     */
    typedef enum OutFileFormat {
        NONE    = 0,
        GMSH    = 1,
        VTK     = 2,
    } OutFileFormat;

    /**
     * Types of reference data
     */
    typedef enum RefType {
        NODE_DATA   = 0,
        CORNER_DATA = 1,
        ELEM_DATA   = 3
    } RefType;

    /**
     * \brief This method returns pointer at existing data, when corresponding
     * output data exists or it creates new one.
     */
    OutputDataBase *output_data_by_field(FieldCommonBase *field, RefType ref_type);

    OutFileFormat   file_format;
    //OutputFormat    *output_format;
    string          *name;              ///< Name of output stream

    /**
     * \brief Vector of pointers at OutputTime
     */
    static std::vector<OutputTime*> output_streams;

    /**
     * \brief Try to find output stream with this name
     */
    static OutputTime *output_stream_by_name(string name);

    /**
     * \brief Does OutputStream with same name and filename exist?
     *
     * When this record is already created, then it returns pointer at
     * corresponding OutputTime. When this record doesn't exist, then
     * it create new OutputTime object and it puts this object to the
     * array of OutputTime pointers
     *
     * \param[in] in_rec  The reference at the input record
     */
    static OutputTime *output_stream(const Input::Record &in_rec);

    /**
     * \brief This method delete all object instances of class OutputTime stored
     * in output_streams vector
     */
    static void destroy_all(void);

    /**
     * \brief This method set current time for registered data array/vector
     */
    void set_data_time(void *data, double time);

    /**
     * \brief The specification of output stream
     *
     * \return This variable defines record for output stream
     */
    static Input::Type::Record input_type;

    /**
     * \brief The specification of output file format
     */
    static Input::Type::AbstractRecord input_format_type;

    int              current_step;      ///< Current step

    /**
     * \brief This method write all registered data to output streams
     */
    static void write_all_data(void);


    /**
     * \brief Generic method for registering output data stored in MultiField
     */
    template<int spacedim, class Value>
    static void register_data(const Input::Record &in_rec,
            const RefType type,
            MultiField<spacedim, Value> *multi_field);


    /**
     * \brief Generic method for registering of output data stored in Field
     *
     * This
     */
    template<int spacedim, class Value>
    static void register_data(const Input::Record &in_rec,
            const RefType ref_type,
            Field<spacedim, Value> *field);

    /**
     * \brief Method for clearing all registered data
     */
    static void clear_data(void);

    /**
     *
     */
    virtual int write_data(void) = 0;

    /**
     *
     */
    virtual int write_head(void) = 0 ;

    /**
     *
     */
    virtual int write_tail(void) = 0;

    /**
     *
     */
    static OutputTime* create_output_stream(const Input::Record &in_rec);

    string *base_filename() { return this->_base_filename; };

private:
    ofstream        *base_file;         ///< Base output stream
    string          *_base_filename;     ///< Name of base output file
    string          *data_filename;     ///< Name of data output file
    ofstream        *data_file;         ///< Data output stream (could be same as base_file)
    Mesh            *mesh;

protected:

    int             rank;               ///< MPI rank of process (is tested in methods)

    // Protected setters for descendant
    void set_mesh(Mesh *_mesh) { mesh = _mesh; };

    void set_base_file(ofstream *_base_file) { this->base_file = _base_file; };

};


template<int spacedim, class Value>
void OutputTime::register_data(const Input::Record &in_rec,
        const RefType type,
        MultiField<spacedim, Value> *multi_field)
{
    for (unsigned long index=0; index < multi_field->size(); index++) {
        OutputTime::register_data(in_rec, type, &multi_field[index]);
    }
}


template<int spacedim, class Value>
void OutputTime::register_data(const Input::Record &in_rec,
        const RefType ref_type,
        Field<spacedim, Value> *field)
{
    string name_ = field->name();
    OutputDataBase *output_data;
    unsigned int item_count = 0, comp_count = 0, node_id;

    // TODO: do not try to find empty string and raise exception

    // Try to find record with output stream (the key is name of data)
    Input::Iterator<string> stream_name_iter = in_rec.find<string>(name_);

    // If record was not found, then throw exception
    if(!stream_name_iter) {
        return;
    }

    // Try to find existing output stream
    OutputTime *output_time = OutputTime::output_stream_by_name(*stream_name_iter);

    /* It's possible now to do output to the file only in the first process */
    if(output_time == NULL || output_time->rank != 0) {
        /* TODO: do something, when support for Parallel VTK is added */
        return;
    }

    Mesh *mesh = (Mesh*)field->mesh();

    if(output_time->get_mesh() == NULL) {
        output_time->set_mesh(mesh);
    }

    ElementFullIter ele = ELEMENT_FULL_ITER(mesh, NULL);
    Node *node;
    int corner_index = 0;
    int node_index = 0;
    int ele_index = 0;

    output_data = output_time->output_data_by_field((FieldCommonBase*)field,
            ref_type);

    switch(ref_type) {
    case NODE_DATA:
        item_count = mesh->n_nodes();
        break;
    case CORNER_DATA:
        // Compute number of all corners
        item_count = 0;
        FOR_ELEMENTS(mesh, ele) {
            item_count += ele->n_nodes();
        }
        break;
    case ELEM_DATA:
        item_count = mesh->n_elements();
        break;
    }

    if(output_data == NULL) {
        /* This is problematic part, because of templates :-( */
        if(typeid(Value) == typeid(FieldValue<1>::Integer)) {
            output_data = (OutputDataBase*)new OutputData<int>(field, item_count, 1);
        } else if(typeid(Value) == typeid(FieldValue<1>::IntVector)) {
            output_data = (OutputDataBase*)new OutputData<int>(field, item_count, 3);
        } else if(typeid(Value) == typeid(FieldValue<1>::Enum)) {
            output_data = (OutputDataBase*)new OutputData<unsigned int>(field, item_count, 1);
        } else if(typeid(Value) == typeid(FieldValue<1>::EnumVector)) {
            output_data = (OutputDataBase*)new OutputData<unsigned int>(field, item_count, 3);
        } else if(typeid(Value) == typeid(FieldValue<1>::Scalar)) {
            output_data = (OutputDataBase*)new OutputData<double>(field, item_count, 1);
        } else if(typeid(Value) == typeid(FieldValue<1>::Vector)) {
            output_data = (OutputDataBase*)new OutputData<double>(field, item_count, 3);
        } else {
            /* TODO: this will not be necessary */
            throw "Try to register unsupported data type.";
        }

        switch(ref_type) {
        case NODE_DATA:
            output_time->node_data.push_back(output_data);
            break;
        case CORNER_DATA:
            output_time->corner_data.push_back(output_data);
            break;
        case ELEM_DATA:
            output_time->elem_data.push_back(output_data);
            break;
        }
    }

    unsigned int *count = new unsigned int[item_count];

    /* Copy data to array */
    switch(ref_type) {
    case NODE_DATA:

        // Initialize arrays
        for(node_id=0; node_id<item_count; node_id++) {
            if(typeid(Value) == typeid(FieldValue<1>::Integer)) {
                (*(OutputData<int>*)output_data)[node_id] = 0;
            } else if(typeid(Value) == typeid(FieldValue<1>::IntVector)) {
                (*(OutputData<int>*)output_data)[node_id] = 0;
            } else if(typeid(Value) == typeid(FieldValue<1>::Enum)) {
                (*(OutputData<unsigned int>*)output_data)[node_id] = 0;
            } else if(typeid(Value) == typeid(FieldValue<1>::EnumVector)) {
                (*(OutputData<unsigned int>*)output_data)[node_id] = 0;
            } else if(typeid(Value) == typeid(FieldValue<1>::Scalar)) {
                (*(OutputData<double>*)output_data)[node_id] = 0;
            } else if(typeid(Value) == typeid(FieldValue<1>::Vector)) {
                (*(OutputData<double>*)output_data)[node_id] = 0;
            }
            count[node_id] = 0;
        }

        /* Copy data to temporary array */
        FOR_ELEMENTS(mesh, ele) {
            FOR_ELEMENT_NODES(ele, node_id) {
                node = ele->node[node_id];
                node_index = mesh->node_vector.index(ele->node[node_id]);
                if(typeid(Value) == typeid(FieldValue<1>::Integer)) {
                    (*(OutputData<int>*)output_data)[node_index] += field->value(node->point(), mesh->element_accessor(ele_index));
                } else if(typeid(Value) == typeid(FieldValue<1>::IntVector)) {
                    (*(OutputData<int>*)output_data)[node_index] += field->value(node->point(), mesh->element_accessor(ele_index));
                } else if(typeid(Value) == typeid(FieldValue<1>::Enum)) {
                    (*(OutputData<unsigned int>*)output_data)[node_index] += field->value(node->point(), mesh->element_accessor(ele_index));
                } else if(typeid(Value) == typeid(FieldValue<1>::EnumVector)) {
                    (*(OutputData<unsigned int>*)output_data)[node_index] += field->value(node->point(), mesh->element_accessor(ele_index));
                } else if(typeid(Value) == typeid(FieldValue<1>::Scalar)) {
                    (*(OutputData<double>*)output_data)[node_index] += field->value(node->point(), mesh->element_accessor(ele_index));
                } else if(typeid(Value) == typeid(FieldValue<1>::Vector)) {
                    (*(OutputData<double>*)output_data)[node_index] += field->value(node->point(), mesh->element_accessor(ele_index));
                }
                count[mesh->node_vector.index(ele->node[node_id])]++;
            }
        }

        // Compute mean values at nodes
        for(node_id=0; node_id<item_count; node_id++) {
            if(typeid(Value) == typeid(FieldValue<1>::Integer)) {
                (*(OutputData<int>*)output_data)[node_id] /= count[node_id];
            } else if(typeid(Value) == typeid(FieldValue<1>::IntVector)) {
                (*(OutputData<int>*)output_data)[node_id] /= count[node_id];
            } else if(typeid(Value) == typeid(FieldValue<1>::Enum)) {
                (*(OutputData<unsigned int>*)output_data)[node_id] /= count[node_id];
            } else if(typeid(Value) == typeid(FieldValue<1>::EnumVector)) {
                (*(OutputData<unsigned int>*)output_data)[node_id] /= count[node_id];
            } else if(typeid(Value) == typeid(FieldValue<1>::Scalar)) {
                (*(OutputData<double>*)output_data)[node_id] /= count[node_id];
            } else if(typeid(Value) == typeid(FieldValue<1>::Vector)) {
                (*(OutputData<double>*)output_data)[node_id] /= count[node_id];
            }
        }

        break;
    case CORNER_DATA:
        FOR_ELEMENTS(mesh, ele) {
            FOR_ELEMENT_NODES(ele, node_id) {
                node = ele->node[node_id];
                if(typeid(Value) == typeid(FieldValue<1>::Integer)) {
                    (*(OutputData<int>*)output_data)[corner_index] = field->value(node->point(), mesh->element_accessor(ele_index));
                } else if(typeid(Value) == typeid(FieldValue<1>::IntVector)) {
                    (*(OutputData<int>*)output_data)[corner_index] = field->value(node->point(), mesh->element_accessor(ele_index));
                } else if(typeid(Value) == typeid(FieldValue<1>::Enum)) {
                    (*(OutputData<unsigned int>*)output_data)[corner_index] = field->value(node->point(), mesh->element_accessor(ele_index));
                } else if(typeid(Value) == typeid(FieldValue<1>::EnumVector)) {
                    (*(OutputData<unsigned int>*)output_data)[corner_index] = field->value(node->point(), mesh->element_accessor(ele_index));
                } else if(typeid(Value) == typeid(FieldValue<1>::Scalar)) {
                    (*(OutputData<double>*)output_data)[corner_index] = field->value(node->point(), mesh->element_accessor(ele_index));
                } else if(typeid(Value) == typeid(FieldValue<1>::Vector)) {
                    (*(OutputData<double>*)output_data)[corner_index] = field->value(node->point(), mesh->element_accessor(ele_index));
                }
                corner_index++;
            }
            ele_index++;
        }
        break;
    case ELEM_DATA:
        FOR_ELEMENTS(mesh, ele) {
            if(typeid(Value) == typeid(FieldValue<1>::Integer)) {
                (*(OutputData<int>*)output_data)[ele_index] = field->value(ele->centre(), mesh->element_accessor(ele_index));
            } else if(typeid(Value) == typeid(FieldValue<1>::IntVector)) {
                (*(OutputData<int>*)output_data)[ele_index] = field->value(ele->centre(), mesh->element_accessor(ele_index));
            } else if(typeid(Value) == typeid(FieldValue<1>::Enum)) {
                (*(OutputData<unsigned int>*)output_data)[ele_index] = field->value(ele->centre(), mesh->element_accessor(ele_index));
            } else if(typeid(Value) == typeid(FieldValue<1>::EnumVector)) {
                (*(OutputData<unsigned int>*)output_data)[ele_index] = field->value(ele->centre(), mesh->element_accessor(ele_index));
            } else if(typeid(Value) == typeid(FieldValue<1>::Scalar)) {
                (*(OutputData<double>*)output_data)[ele_index] = field->value(ele->centre(), mesh->element_accessor(ele_index));
            } else if(typeid(Value) == typeid(FieldValue<1>::Vector)) {
                (*(OutputData<double>*)output_data)[ele_index] = field->value(ele->centre(), mesh->element_accessor(ele_index));
            }
            ele_index++;
        }
        break;
    }

    /* Set the last time */
    if(output_time->time < field->time()) {
        output_time->time = field->time();
    }
}


#endif
