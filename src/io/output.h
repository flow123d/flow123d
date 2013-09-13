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
#include <mpi.h>
#include <boost/any.hpp>

#include "system/xio.h"
#include "mesh/mesh.h"

#include "fields/field_base.hh"
#include "input/accessors.hh"

class OutputFormat;

/**
 * \brief This method is used for stored data that are copied from field.
 */
class OutputData {
public:
    FieldCommonBase *field;
    std::vector<boost::any> data;
    OutputData(FieldCommonBase *field);
    ~OutputData();
};


/**
 * \brief The class for outputing data during time.
 *
 * This class is descendant of Output class. This class is used for outputing
 * data varying in time. Own output to specific file formats is done at other
 * places to. See output_vtk.cc and output_msh.cc.
 */
class OutputTime {
protected:

    int             rank;               ///< MPI rank of process (is tested in methods)

    // Protected setters for descendant
    void set_mesh(Mesh *_mesh) { mesh = _mesh; };

    void set_base_file(ofstream *_base_file) { base_file = _base_file; };

    void set_base_filename(string *_base_filename) { base_filename = _base_filename; };

    OutputTime() {};

private:
    ofstream        *base_file;         ///< Base output stream
    string          *base_filename;     ///< Name of base output file
    string          *data_filename;     ///< Name of data output file
    ofstream        *data_file;         ///< Data output stream (could be same as base_file)
    Mesh            *mesh;

public:
    vector<OutputData*>    node_data;
    vector<OutputData*>    corner_data;
    vector<OutputData*>    elem_data;

    double          time;               ///< The newest time of registered data

    double          write_time;         ///< The last time, when data was wrote to this stream

    ofstream& get_base_file(void) { return *base_file; };

    string& get_base_filename(void) { return *base_filename; };

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

    OutFileFormat   file_format;
    OutputFormat    *output_format;
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
    void register_data(const Input::Record &in_rec,
            const RefType type,
            Field<spacedim, Value> *field);

    /**
     * \brief Method for clearing all registered data
     */
    static void clear_data(void);

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
        const RefType type,
        Field<spacedim, Value> *field)
{
    string name_ = field->name();
    OutputData *output_data;

    // Try to find record with output stream (the key is name of data)
    Input::Iterator<string> stream_name_iter = in_rec.find<string>(name_);

    // If record was not found, then exit
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

    Mesh *mesh = output_time->get_mesh();
    ElementFullIter ele = ELEMENT_FULL_ITER(mesh, NULL);

    /* Copy data to vector */
    switch(type) {
    case NODE_DATA:
        output_data = new OutputData((FieldCommonBase*)field);
        output_time->node_data.push_back(output_data);
        break;
    case CORNER_DATA:
        output_data = new OutputData((FieldCommonBase*)field);
        output_time->corner_data.push_back(output_data);
        break;
    case ELEM_DATA:
        output_data = new OutputData((FieldCommonBase*)field);
        FOR_ELEMENTS(mesh, ele) {
            output_data->data.push_back(field->value(ele->centre(), ele->element_accessor()));
        }
        output_time->elem_data.push_back(output_data);
        break;
    }

    /* Set the last time */
    if(output_time->time < field->time()) {
        output_time->time = field->time();
    }
}


/**
 * \brief The class used as parent class of file format classes
 */
class OutputFormat {
public:
	OutputFormat() {}
    virtual ~OutputFormat() {}
	virtual int write_data(void) { return 0; }
	virtual int write_head(void) { return 0; }
	virtual int write_tail(void) { return 0; }

	static Input::Type::AbstractRecord input_type;
};


#endif
