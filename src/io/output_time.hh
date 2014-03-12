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

//#include <mpi.h>

#include "mesh/mesh.h"
#include "input/accessors.hh"



class OutputDataBase;
class FieldCommonBase; // in fact not necessary, output_data_by_field() can use directly name as parameter
template <int spacedim, class Value>
class Field;
template <int spacedim, class Value>
class MultiField;


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
     * \brief Try to find output stream from a key in record.
     */
    static OutputTime *output_stream_by_key_name(const Input::Record &in_rec, const string key_name);

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
            MultiField<spacedim, Value> &multi_field);


    /**
     * \brief Generic method for registering of output data stored in Field
     *
     * This
     */
    template<int spacedim, class Value>
    static void register_data(const Input::Record &in_rec,
            const RefType ref_type,
            Field<spacedim, Value> &field);

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

    template<int spacedim, class Value>
    static void compute_field_data(const RefType ref_type, Field<spacedim, Value> *field, OutputTime *output_time);


private:
    ofstream        *base_file;         ///< Base output stream
    string          *_base_filename;     ///< Name of base output file
    string          *data_filename;     ///< Name of data output file
    ofstream        *data_file;         ///< Data output stream (could be same as base_file)
    Mesh      *mesh;

protected:

    int             rank;               ///< MPI rank of process (is tested in methods)

    // Protected setters for descendant
    void set_mesh(Mesh *_mesh) { mesh = _mesh; };

    void set_base_file(ofstream *_base_file) { this->base_file = _base_file; };

};




#endif /* OUTPUT_TIME_HH_ */
