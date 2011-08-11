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
 */

#ifndef OUTPUT_H
#define OUTPUT_H

#include "mesh/mesh.h"

#include <vector>
#include <string>
#include <fstream>

// TODO: move enums to classes

// FILE formats
typedef enum {
    GMSH_MSH_ASCII = 1,
    GMSH_MSH_BIN = 2,
    VTK_SERIAL_ASCII = 3,
    VTK_PARALLEL_ASCII = 4
} OutFileFormat;

// Types of data, that could be written to output file
typedef enum {
    OUT_VECTOR_INT_SCA,
    OUT_VECTOR_INT_VEC,
    OUT_VECTOR_FLOAT_SCA,
    OUT_VECTOR_FLOAT_VEC,
    OUT_VECTOR_DOUBLE_SCA,
    OUT_VECTOR_DOUBLE_VEC,
    OUT_ARRAY_INT_SCA,
    OUT_ARRAY_FLOAT_SCA,
    OUT_ARRAY_DOUBLE_SCA
} OutDataType;

/**
 * Class of output data storing reference on data.
 *
 * This class is referenced in Output and OuputTime class. The object contains
 * only reference on data. The referenced data could not be freed or rewrite
 * until data are written to the output file.
 */
class OutputData {
private:
public:
    string          *name;      ///< String with name of data
    string          *units;     ///< String with units
    void            *data;      ///< Pointer at own data
    OutDataType     type;       ///< Type values in vector
    int             comp_num;   ///< Number of components in vector
    int             num;        ///< Number of values in vector/array

    OutputData() {};            ///< Un-named constructor can't be called directly

    string* getName(void) { return name; };
    string* getUnits(void) { return units; };
    int getCompNum(void) { return comp_num; };
    int getValueNum(void) { return num; };

    /**
     * \brief Constructor for OutputData storing names of output data and their
     * units.
     */
    OutputData(std::string name, std::string unit, int *data, unsigned int size);

    /**
     * \brief Constructor for OutputData storing names of output data and their
     * units.
     */
    OutputData(std::string name, std::string unit, float *data, unsigned int size);

    /**
     * \brief Constructor for OutputData storing names of output data and their
     * units.
     */
    OutputData(std::string name, std::string unit, double *data, unsigned int size);

    /**
     * \brief Constructor for OutputData storing names of output data and their
     * units.
     */
    OutputData(std::string name, std::string unit, std::vector<int> &data);

    /**
     * \brief Constructor for OutputData storing names of output data and their
     * units.
     */
    OutputData(std::string name, std::string unit, std::vector< vector<int> > &data);

    /**
     * \brief Constructor for OutputData storing names of output data and their
     * units.
     */
    OutputData(std::string name, std::string unit, std::vector<float> &data);

    /**
     * \brief Constructor for OutputData storing names of output data and their
     * units.
     */
    OutputData(std::string name, std::string unit, std::vector< vector<float> > &data);

    /**
     * \brief Constructor for OutputData storing names of output data and their
     * units.
     */
    OutputData(std::string name, std::string unit, std::vector<double> &data);

    /**
     * \brief Constructor for OutputData storing names of output data and their
     * units.
     */
    OutputData(std::string name, std::string unit, std::vector< vector<double> > &data);

    /**
     * \brief Destructor for OutputData
     */
    ~OutputData();
};

/**
 * Definition of output data vector
 */
typedef std::vector<OutputData> OutputDataVec;

/**
 * Class of output. This class is used for output data to the output files.
 *
 * This class is used for writing static data to the output files. This class
 * itself doesn't include any method for writing data to the output file.
 * The methods for writing data to the output files are in places. Look at
 * output_vtk.cc and output_msh.cc that implement output for specific file
 * formats.
 */
class Output {
public:
    /*
     * \brief Temporary definition for storing data (C++ vector of double scalars)
     */
    typedef std::vector<double> ScalarFloatVector;

    /**
     * \brief Temporary definition for storing data (C++ vector of double vectors)
     */
    typedef std::vector< vector<double> > VectorFloatVector;

    /**
     * \brief Temporary structure for storing data (double scalars)
     */
    typedef struct OutScalar {
        ScalarFloatVector   *scalars;
        string              name;
        string              unit;
    } _OutScalar;

    /**
     * \brief Temporary structure for storing data (double vectors)
     */
    typedef struct OutVector {
        VectorFloatVector   *vectors;
        string              name;
        string              unit;
    } _OutVector;

    /**
     * \brief Temporary vectors of structures for storing data (double scalars)
     */
    typedef std::vector<OutScalar> OutScalarsVector;

    /**
     * \brief Temporary vectors of structures for storing data (double vectors)
     */
    typedef std::vector<OutVector> OutVectorsVector;

    /**
     * \brief Constructor of the Output object
     *
     * \param[in] *_mesh    The pointer at Mesh
     * \param[in] *fname    The name of the output file
     */
    Output(Mesh *mesh, string filename);

    /**
     * Convert string output format names into enum values.
     */
    OutFileFormat parse_output_format(char* format_name);

    /**
     * \brief Destructor of the Output object.
     */
    ~Output();

    /**
     * \brief This function gets data from mesh and save them in Output. It is temporary solution.
     */
    void get_data_from_mesh(void);

    /**
     * \brief This function free data from Mesh. It is temporary solution.
     */
    void free_data_from_mesh(void);

    /**
     * \brief Register array of data on nodes.
     *
     * This function will add reference on the array of data to the Output object.
     * Own data will be written to the file, when write_data() method will be called.
     *
     * \param[in]   name    The name of data
     * \param[in]   unit    The name of units
     * \param[in]   data    The pointer on array of data
     * \param[in]   size    The number of values in array
     *
     * \return This function return 1, when the length of data vector is the same as
     * number of nodes in mesh. When the the number is different, then this
     * function returns 0.
     */
    template <typename _Data>
    int register_node_data(std::string name, std::string unit, _Data *data, uint size);

    /**
     * \brief Register array of data on elements.
     *
     * This function will add reference on this array of data to the Output object.
     * Own data will be written to the file, when write_data() method will be called.
     *
     * \param[in]   name    The name of data
     * \param[in]   unit    The name of units
     * \param[in]   data    The pointer on array of data
     * \param[in]   size    The number of values in array
     *
     * \return This function return 1, when the length of data vector is the same as
     * number of elements in mesh. When the the number is different, then this
     * function returns 0.
     */
    template <typename _Data>
    int register_elem_data(std::string name, std::string unit, _Data *data, unsigned int size);

    /**
     * \brief Register data on nodes.
     *
     * This function will add reference on this data to the Output object. Own data
     * will be written to the file, when write_data() method will be called.
     *
     * \param[in]   name    The name of data
     * \param[in]   unit    The name of units
     * \param[in]   data    The reference on vector of data
     *
     * \return This function return 1, when the length of data vector is the same as
     * number of nodes in mesh. When the the number is different, then this function
     * returns 0.
     */
    template <typename _Data>
    int register_node_data(std::string name, std::string unit, std::vector<_Data> &data);

    /**
     * \brief Register data on elements.
     *
     * This function will add reference on the data to the Output object. Own data
     * will be written to the file, when write_data() method will be called.
     *
     * \param[in]   name    The name of data
     * \param[in]   unit    The name of units
     * \param[in]   data    The reference on vector of data
     *
     * \return This function return 1, when the length of data vector is the same as
     * number of elements in mesh. When the the number is different, then this
     * function returns 0.
     */
    template <typename _Data>
    int register_elem_data(std::string name, std::string unit, std::vector<_Data> &data);

    /**
     * \brief This function call pointer at _write_data(Output). It writes
     * registered data to specified file format.
     *
     * \return This function return result of pointer at output function.
     */
    int write_data(void);

    // TODO: Move public getters to protected section
    std::vector<OutputData> *get_node_data(void) { return node_data; };
    std::vector<OutputData> *get_elem_data(void) { return elem_data; };
    ofstream& get_base_file(void) { return *base_file; };
    string& get_base_filename(void) { return *base_filename; };
    ofstream& get_data_file(void) { return *data_file; };
    string& get_data_filename(void) { return *data_filename; };
    Mesh *get_mesh(void) { return mesh; };
    char get_format_type(void) { return format_type; };

    // TODO: Move public setters to protected section
    void set_data_file(ofstream *_data_file) { data_file = _data_file; };

protected:
    // Protected setters for descendant
    void set_mesh(Mesh *_mesh) { mesh = _mesh; };
    void set_base_file(ofstream *_base_file) { base_file = _base_file; };
    void set_base_filename(string *_base_filename) { base_filename = _base_filename; };
    void set_format_type(OutFileFormat _format_type) { format_type = _format_type; };
    void set_node_data(std::vector<OutputData> *_node_data) { node_data = _node_data; };
    void set_elem_data(std::vector<OutputData> *_elem_data) { elem_data = _elem_data; };

    Output() { node_scalar = NULL; element_scalar = NULL; element_vector = NULL; };

private:
    struct OutScalar *node_scalar;      // Temporary solution
    struct OutScalar *element_scalar;   // Temporary solution
    struct OutVector *element_vector;   // Temporary solution

    ofstream        *base_file;         ///< Base output stream
    string          *base_filename;     ///< Name of base output file
    string          *data_filename;     ///< Name of data output file
    ofstream        *data_file;         ///< Data output stream (could be same as base_file)
    OutFileFormat             format_type;        ///< Type of output
    Mesh            *mesh;
    OutputDataVec   *node_data;         ///< List of data on nodes
    OutputDataVec   *elem_data;         ///< List of data on elements

    // TODO: tohle by se melo resit jinak. napr: Output by melo obsahovat
    // abstraktni virtualni metodu write_data
    // pak bychom meli potomky tridy Output pro jednotlive formaty vystupu,
    // a ty by definovali vystup do patricneho formatu

    // Internal API for file formats
    int (*_write_data)(Output *output);
};

template <typename _Data>
int Output::register_node_data(std::string name,
        std::string unit,
        _Data *data,
        uint size)
{
	int rank=0;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	/* It's possible now to do output to the file only in the first process */
	if(rank!=0) {
		/* TODO: do something, when support for Parallel VTK is added */
		return 0;
	}

    if(mesh->node_vector.size() == size) {
        int found = 0;

        OutputData *out_data = new OutputData(name, unit, data, size);
        node_data->push_back(*out_data);

        return 1;
    } else {
        xprintf(Err, "Number of values: %d is not equal to number of nodes: %d\n", size, mesh->node_vector.size());
        return 0;
    }
}

template <typename _Data>
int Output::register_elem_data(std::string name,
        std::string unit,
        _Data *data,
        uint size)
{
	int rank=0;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	/* It's possible now to do output to the file only in the first process */
	if(rank!=0) {
		/* TODO: do something, when support for Parallel VTK is added */
		return 0;
	}

    if(mesh->element.size() == size) {
        int found = 0;

        OutputData *out_data = new OutputData(name, unit, data, size);
        elem_data->push_back(*out_data);

        return 1;
    } else {
        xprintf(Err, "Number of values: %d is not equal to number of elements: %d\n", size, mesh->element.size());
        return 0;
    }
}

template <typename _Data>
int Output::register_node_data(std::string name,
        std::string unit,
        std::vector<_Data> &data)
{
	int rank=0;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	/* It's possible now to do output to the file only in the first process */
	if(rank!=0) {
		/* TODO: do something, when support for Parallel VTK is added */
		return 0;
	}

    if(mesh->node_vector.size() == data.size()) {
        OutputData *out_data = new OutputData(name, unit, data);
        node_data->push_back(*out_data);
        return 1;
    } else {
        xprintf(Err, "Number of values: %d is not equal to number of nodes: %d\n", data.size(), mesh->node_vector.size());
        return 0;
    }
}

template <typename _Data>
int Output::register_elem_data(std::string name,
        std::string unit,
        std::vector<_Data> &data)
{
	int rank=0;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	/* It's possible now to do output to the file only in the first process */
	if(rank!=0) {
		/* TODO: do something, when support for Parallel VTK is added */
		return 0;
	}

    if(mesh->element.size() == data.size()) {
        OutputData *out_data = new OutputData(name, unit, data);
        elem_data->push_back(*out_data);
        return 1;
    } else {
        xprintf(Err, "Number of values: %d is not equal to number of elements: %d\n", data.size(), mesh->element.size());
        return 0;
    }
}


/**
 * \brief The class for outputing data during time.
 *
 * This class is descendant of Output class. This class is used for outputing
 * data varying in time. Own output to specific file formats is done at other
 * places to. See output_vtk.cc and output_msh.cc.
 */
class OutputTime : public Output {
public:
    /**
     * \brief Constructor of OutputTime object. It opens base file for writing.
     *
     * \param[in]   *_mesh  The pointer at mesh object.
     * \param[in]   fname   The name of output file
     */
    OutputTime(Mesh *mesh, string filename);

    /**
     * \brief Destructor of OutputTime. It doesn't do anything, because all
     * necessary destructors will be called in destructor of Output
     */
    ~OutputTime();

    /**
     * \brief This function register data on nodes.
     *
     * This function will add reference on this array of data to the Output object.
     * It is possible to call this function only once, when data are at the same
     * address during time. It is possible to call this function for each step, when
     * data are not at the same address, but name of the data has to be same.
     * Own data will be written to the file, when write_data() method will be called.
     *
     * \param[in] name  The name of data
     * \param[in] unit  The units of data
     * \param[in] *data The pointer at data (array of int, float or double)
     * \param[in] size  The size of array (number of values)
     *
     * \return This function returns 1, when data were registered. This function
     * returns 0, when it wasn't able to register data (number of values isn't
     * same as number of nodes).
     */
    template <typename _Data>
    int register_node_data(std::string name, std::string unit, _Data *data, unsigned int size);

    /**
     * \brief This function register data on elements.
     *
     * This function will add reference on this array of data to the Output object.
     * It is possible to call this function only once, when data are at the same
     * address during time. it is possible to call this function for each step, when
     * data are not at the same address, but name of the data has to be same.
     * Own data will be written to the file, when write_data() method will be called.
     *
     * \param[in] name  The name of data
     * \param[in] unit  The units of data
     * \param[in] *data The pointer at data (array of int, float or double)
     * \param[in] size  The size of array (number of values)
     *
     * \return This function returns 1, when data were registered. This function
     * returns 0, when it wasn't able to register data (number of values isn't
     * same as number of elements).
     */
    template <typename _Data>
    int register_elem_data(std::string name, std::string unit, _Data *data, unsigned int size);

    /**
     * \brief This function register data on nodes.
     *
     * This function will add reference on this array of data to the Output object.
     * It is possible to call this function only once, when data are at the same
     * address during time. it is possible to call this function for each step, when
     * data are not at the same address, but name of the data has to be same.
     * Own data will be written to the file, when write_data() method will be called.
     *
     * \param[in] name  The name of data
     * \param[in] unit  The units of data
     * \param[in] *data The pointer at data (array of int, float or double)
     * \param[in] size  The size of array (number of values)
     *
     * \return This function returns 1, when data were registered. This function
     * returns 0, when it wasn't able to register data (number of values isn't
     * same as number of nodes).
     */
    template <typename _Data>
    int register_node_data(std::string name, std::string unit, std::vector<_Data> &data);

    /**
     * \brief Register vector of data on elements.
     *
     * This function will add reference on the data to the Output object. Own
     * data will be written to the file, when write_data() method will be called.
     * When the data has been already registered, then pointer at data will be
     * updated. Otherwise, new data will be registered.
     *
     * \param[in] name  The name of data
     * \param[in] unit  The unit of data
     * \param[in] &data The reference on vector (int, float, double)
     *
     * \return This function returns 1, when data were successfully registered.
     * This function returns 0, when number of elements and items of vector is
     * not the same.
     */
    template <typename _Data>
    int register_elem_data(std::string name, std::string unit, std::vector<_Data> &data);

    /**
     * \brief This function call pointer at appropriate pointer at function,
     * that write data to specific file format.
     *
     * \param[in] time  The output will be done for this time
     *
     * \return This function returns result of method _write_data().
     */
    int write_data(double time);
private:
    int              current_step;      ///< Current step
    // Internal API for file formats
    int (*_write_data)(OutputTime *output, double time, int step);
    int (*_write_head)(OutputTime *output);
    int (*_write_tail)(OutputTime *output);
};

template <typename _Data>
int OutputTime::register_node_data(std::string name,
        std::string unit,
        _Data *data,
        uint size)
{
	int rank=0;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	int found = 0;

	/* It's possible now to do output to the file only in the first process */
	if(rank!=0) {
		/* TODO: do something, when support for Parallel VTK is added */
		return 0;
	}

    std::vector<OutputData> *node_data = get_node_data();
    Mesh *mesh = get_mesh();

    ASSERT(mesh->node_vector.size() == size,
            "mesh->node_vector.size(): %d != size: %d",
            mesh->node_vector.size(),
            size);

    for(std::vector<OutputData>::iterator od_iter = node_data->begin();
            od_iter != node_data->end();
            od_iter++)
    {
        if(*od_iter->name == name) {
            od_iter->data = (void*)data;
            found = 1;
            break;
        }
    }

    if(found == 0) {
        OutputData *out_data = new OutputData(name, unit, data, size);
        node_data->push_back(*out_data);
    }

    return 1;

}

template <typename _Data>
int OutputTime::register_elem_data(std::string name,
        std::string unit,
        _Data *data,
        unsigned int size)
{
	int rank=0;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    int found = 0;

	/* It's possible now to do output to the file only in the first process */
	if(rank!=0) {
		/* TODO: do something, when support for Parallel VTK is added */
		return 0;
	}

    std::vector<OutputData> *elem_data = get_elem_data();
    Mesh *mesh = get_mesh();

    ASSERT(mesh->element.size() == size,
            "mesh->element.size(): %d != size: %d",
            mesh->element.size(),
            size);

    for(std::vector<OutputData>::iterator od_iter = elem_data->begin();
            od_iter != elem_data->end();
            od_iter++)
    {
        if(*od_iter->name == name) {
            od_iter->data = (void*)data;
            found = 1;
            break;
        }
    }

    if(found == 0) {
        OutputData *out_data = new OutputData(name, unit, data, size);
        elem_data->push_back(*out_data);
    }

    return 1;
}

template <typename _Data>
int OutputTime::register_node_data(std::string name,
        std::string unit,
        std::vector<_Data> &data)
{
	int rank=0;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    int found = 0;

	/* It's possible now to do output to the file only in the first process */
	if(rank!=0) {
		/* TODO: do something, when support for Parallel VTK is added */
		return 0;
	}

    std::vector<OutputData> *node_data = get_node_data();
    Mesh *mesh = get_mesh();

    ASSERT(mesh->node_vector.size() == data.size(),
            "mesh->node_vector.size(): %d != size: %d",
            mesh->node_vector.size(),
            data.size());

    for(std::vector<OutputData>::iterator od_iter = node_data->begin();
            od_iter != node_data->end();
            od_iter++)
    {
        if(*od_iter->name == name) {
            od_iter->data = (void*)&data;
            found = 1;
            break;
        }
    }

    if(found == 0) {
        OutputData *out_data = new OutputData(name, unit, data);
        node_data->push_back(*out_data);
    }

    return 1;

}

template <typename _Data>
int OutputTime::register_elem_data(std::string name,
        std::string unit,
        std::vector<_Data> &data)
{
	int rank=0;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    int found = 0;

	/* It's possible now to do output to the file only in the first process */
	if(rank!=0) {
		/* TODO: do something, when support for Parallel VTK is added */
		return 0;
	}

    std::vector<OutputData> *elem_data = get_elem_data();
    Mesh *mesh = get_mesh();

    ASSERT(mesh->element.size() == data.size(),
            "mesh->element.size(): %d != size: %d",
            mesh->element.size(),
            data.size());

    for(std::vector<OutputData>::iterator od_iter = elem_data->begin();
            od_iter != elem_data->end();
            od_iter++)
    {
        if(*od_iter->name == name) {
            od_iter->data = (void*)&data;
            found = 1;
            break;
        }
    }

    if(found == 0) {
        OutputData *out_data = new OutputData(name, unit, data);
        elem_data->push_back(*out_data);
    }
    return 1;

}

#endif
