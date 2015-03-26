/*
 * output_time.hh
 *
 *  Created on: Mar 12, 2014
 *      Author: jb
 */

#ifndef OUTPUT_TIME_HH_
#define OUTPUT_TIME_HH_

#include <vector>
#include <string>
#include <fstream>

#include "mesh/mesh.h"
#include "input/accessors.hh"


class OutputDataBase;
class FieldCommon; // in fact not necessary, output_data_by_field() can use directly name as parameter
template <int spacedim, class Value>
class Field;
template <int spacedim, class Value>
class MultiField;
class TimeGovernor;

/**
 * \brief The class for outputting data during time.
 *
 * This class is descendant of Output class. This class is used for outputting
 * data varying in time. Own output to specific file formats is done at other
 * places to. See output_vtk.cc and output_msh.cc.
 */
class OutputTime {

public:
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
     * \brief The specification of output stream
     *
     * \return This variable defines record for output stream
     */
    static Input::Type::Record input_type;

    /**
     * \brief The specification of output file format
     */
    static Input::Type::AbstractRecord input_format_type;

    /**
     * Types of reference data
     */
    enum DiscreteSpace {
        NODE_DATA   = 0,
        CORNER_DATA = 1,
        ELEM_DATA   = 3
    };

    /**
     * \brief This method delete all object instances of class OutputTime stored
     * in output_streams vector
     */
    static void destroy_all(void);

    /**
     * \brief This method tries to create new instance of OutputTime according
     * record in configuration file.
     */
    static std::shared_ptr<OutputTime> create_output_stream(const Input::Record &in_rec);

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
     * \brief Registers names of output fields that can be written using this stream.
     * @param in_array Array of admissible fields (array of selections).
     */
    void add_admissible_field_names(const Input::Array &in_array);

    /**
     * \brief Clear data for output computed by method @p compute_field_data.
     */
    void clear_data(void);

    /**
     *  Add time marks matching given @p tg.output_mark_type as well as general output type
     *  TimeMarks::type_output(). The time marks denotes times when output should be performed according
     *  to the input record of the output stream, namely keys: time_step, time_list, and include_input_times.
     */
    void mark_output_times(const TimeGovernor &tg);

    /**
     * Declaration of new exception info used in following exception
     */
    TYPEDEF_ERR_INFO(EI_FieldName, std::string);

    /**
     * Declaration of exception
     */
    DECLARE_EXCEPTION(ExcOutputVariableVector, << "Can not output field " << EI_FieldName::qval
            << " returning variable size vectors. Try convert to MultiField.\n");

    /**
     * Record for current output stream
     */
    Input::AbstractRecord format_record_;

protected:

    /**
     * Interpolate given @p field into output discrete @p space and store the values
     * into private storage for postponed output.
     */
    template<int spacedim, class Value>
    void compute_field_data(DiscreteSpace space, Field<spacedim, Value> &field);

    /**
     * \brief This method returns pointer at existing data, when corresponding
     * output data exists or it creates new one.
     */
    OutputDataBase *output_data_by_field_name(const string &field_name, DiscreteSpace ref_type);

    /**
     * \brief This method set current time for registered data array/vector
     */
    void set_data_time(void *data, double time);

    /**
     * \brief Virtual method for writing header of output file
     */
    virtual int write_head(void) = 0 ;

    /**
     * \brief Virtual method for writing data to output file
     */
    virtual int write_data(void) = 0;

    /**
     * \brief Virtual method for writing tail of output file
     */
    virtual int write_tail(void) = 0;

    /**
     * \brief Vector of pointers at OutputTime
     */
    static std::vector<OutputTime*> output_streams;

    /**
     * Cached MPI rank of process (is tested in methods)
     */
    int rank;

    /**
     * Vector of registered output data related to nodes
     */
    vector<OutputDataBase*> node_data;

    /**
     * Vector of registered output data related to corner of elements
     */
    vector<OutputDataBase*> corner_data;

    /**
     * Vector of registered output data related to elements
     */
    vector<OutputDataBase*> elem_data;

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
     * Map of names of output fields. True means that field will be saved.
     */
    map<string, bool> output_names;

    /**
     * Record for current output stream
     */
    Input::Record input_record_;

    /**
     * Base output stream
     */
    ofstream *_base_file;

    /**
     * Name of base output file
     */
    string _base_filename;

    /**
     * Data output stream (could be same as base_file)
     */
    ofstream *_data_file;

    /**
     * Cached pointer at mesh used by this output stream
     */
    Mesh *_mesh;
};


#endif /* OUTPUT_TIME_HH_ */
