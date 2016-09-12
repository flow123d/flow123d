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
 * @file    output_time.hh
 * @brief   
 */

#ifndef OUTPUT_TIME_HH_
#define OUTPUT_TIME_HH_

#include <vector>
#include <string>
#include <fstream>
#include "input/accessors.hh"

class FilePath;
class Observe;
class OutputDataBase;
class Mesh;
class FieldCommon; // in fact not necessary, output_data_by_field() can use directly name as parameter
template <int spacedim, class Value>
class Field;
class FieldSet;
template <int spacedim, class Value>
class MultiField;
class TimeGovernor;
class OutputMesh;
class OutputMeshDiscontinuous;


/**
 * \brief The class for outputting data during time.
 *
 * This class is descendant of Output class. This class is used for outputting
 * data varying in time. Own output to specific file formats is done at other
 * places to. See output_vtk.cc and output_msh.cc.
 */
class OutputTime {

public:
    /// Default constructor. Only for testing.
    OutputTime();


    /**
     * \brief Constructor of OutputTime object. It opens base file for writing.
     *
     * \param[in] equation_name The name of equation, used for forming output file name.
     * \param[in] in_rec The reference on the input record
     */
    void init_from_input(const std::string &equation_name, Mesh &mesh, const Input::Record &in_rec);

    /**
     * \brief Destructor of OutputTime. It doesn't do anything, because all
     * necessary destructors will be called in destructor of Output
     */
    virtual ~OutputTime();

    /**
     * Return the input array for the output time set of the output stream.
     */
    Input::Iterator<Input::Array> get_time_set_array();

    /**
     * \brief The specification of output stream
     *
     * \return This variable defines record for output stream
     */
    static const Input::Type::Record & get_input_type();

    /**
     * \brief The specification of output file format
     */
    static Input::Type::Abstract & get_input_format_type();

    /**
     * Types of reference data
     */
    static const unsigned int N_DISCRETE_SPACES = 3;
    enum DiscreteSpace {
        NODE_DATA   = 0,
        CORNER_DATA = 1,
        ELEM_DATA   = 2
    };

    /**
     * \brief This method delete all object instances of class OutputTime stored
     * in output_streams vector
     */
    //static void destroy_all(void);

    /**
     * \brief This method tries to create new instance of OutputTime according
     * record in configuration file.
     */
    static std::shared_ptr<OutputTime> create_output_stream(const std::string &equation_name, Mesh &mesh, const Input::Record &in_rec);
    
    /**
     * Create the output mesh from the given computational mesh. The field set passed in is used
     * to select the field used for adaptivity of the output mesh.
     */
    void make_output_mesh(FieldSet &output_fields);
    
    /**
     * \brief Generic method for registering output data stored in MultiField
     *
     * @param ref_type    Type of output (element, node, corner data).
     * @param multi_field The actual field for output.
     */
    template<int spacedim, class Value>
    void register_data(const DiscreteSpace type,
            MultiField<spacedim, Value> &multi_field);

    /**
     * \brief Generic method for registering of output data stored in Field
     *
     * @param ref_type  Type of output (element, node, corner data).
     * @param field     The actual field for output.
     */
    template<int spacedim, class Value>
    void register_data(const DiscreteSpace ref_type,
            Field<spacedim, Value> &field);

    /**
     * Write all data registered as a new time frame.
     */
    void write_time_frame();

    /**
     * Getter of the oubserve object.
     */
    std::shared_ptr<Observe> observe();

    /**
     * \brief Clear data for output computed by method @p compute_field_data.
     */
    void clear_data(void);

    /**
     * Declaration of new exception info used in following exception
     */
    TYPEDEF_ERR_INFO(EI_FieldName, std::string);

    /**
     * Declaration of exception
     */
    DECLARE_EXCEPTION(ExcOutputVariableVector, << "Can not output field " << EI_FieldName::qval
            << " returning variable size vectors. Try convert to MultiField.\n");


protected:
    
    void compute_discontinuous_output_mesh();
    
    /**
     * Interpolate given @p field into output discrete @p space and store the values
     * into private storage for postponed output.
     */
    template<int spacedim, class Value>
    void compute_field_data(DiscreteSpace type, Field<spacedim, Value> &field);

    /**
     * Change main filename to have prescribed extension.
     */
    void fix_main_file_extension(std::string extension);


    /**
     * \brief Virtual method for writing data to output file
     */
    virtual int write_data(void) = 0;

    /**
     * Cached MPI rank of process (is tested in methods)
     */
    int rank;

    /**
     * Map field name to its OutputData object.
     */
    typedef std::shared_ptr<OutputDataBase> OutputDataPtr;
    typedef std::vector< OutputDataPtr > OutputDataFieldVec;

    /**
     * Registered output data. Single map for every value of DiscreteSpace
     * corresponding to nodes, elements and corners.
     */
    OutputDataFieldVec  output_data_vec_[N_DISCRETE_SPACES];

    /**
     * Current step
     */
    int current_step;

    /**
     * The newest time of registered data
     */
    double time;

    /**
     * The last time, when data was wrote to this stream
     */
    double write_time;

    /**
     * Maps names of output fields required by user to their indices in
     * output_data_vec_.
     */
    typedef unsigned int DiscreteSpaceFlags;

    /**
     * Record for current output stream
     */
    Input::Record input_record_;

    /**
     * Base output stream
     */
    ofstream _base_file;

    /**
     * Name of base output file
     */
    FilePath _base_filename;

    /**
     * Name of the equation owning the output stream. Usually the balance equation.
     * Used for forming default output file name and the name of observe output file.
     */
    std::string equation_name_;

    /**
     * Cached pointer at mesh used by this output stream
     */
    Mesh *_mesh;
    
    /// Output mesh.
    std::shared_ptr<OutputMesh> output_mesh_;
    /// Discontinuous (non-conforming) mesh. Used for CORNER_DATA.
    std::shared_ptr<OutputMeshDiscontinuous> output_mesh_discont_;
    
    std::shared_ptr<Observe> observe_;

    /// Auxliary flag for refinement enabling, due to gmsh format.
    bool enable_refinement_;
};


#endif /* OUTPUT_TIME_HH_ */
