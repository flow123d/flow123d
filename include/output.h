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
 * @file
 * @brief   Header: The functions for all outputs.
 *
 */

#ifndef OUTPUT_H
#define OUTPUT_H

#include <vector>
#include <string>
#include <fstream>

//#include "system.hh"
#include "transport.h"  // TODO remove

/// External types
struct Problem;         // TODO remove
class ConvectionTransport;  // TODO remove
class Mesh;

//TODO: v C++ by mely byt konstanty uvnitr definice trid, nejlepe jako enum napr:
// enum OutputDataFormat { GMSH_ASCII, GMSH_BIN, ..}
// a uvazit zda to ma byt public, nebo private (nebo protected)


// FILE formats
#define POS_ASCII           1   //TODO GMSH_ASCII
#define POS_BIN             2   //TODO
#define VTK_SERIAL_ASCII    3
#define VTK_PARALLEL_ASCII  4

// VTK Element types
#define VTK_VERTEX          1
#define VTK_POLY_VERTEX     2
#define VTK_LINE            3
#define VTK_POLY_LINE       4
#define VTK_TRIANGLE        5
#define VTK_TRIANGLE_STRIP  6
#define VTK_POLYGON         7
#define VTK_PIXEL           8
#define VTK_QUAD            9
#define VTK_TETRA           10
#define VTK_VOXEL           11
#define VTK_HEXAHEDRON      12
#define VTK_WEDGE           13
#define VTK_PYRAMID         14

#define VTK_QUADRIC_EDGE        21
#define VTK_QUADRIC_TRIANGLE    22
#define VTK_QUADRIC_QUAD        23
#define VTK_QUADRIC_TETRA       24
#define VTK_QUADRIC_HEXAHEDRON  25

// VTK Element size (number of nodes)
#define VTK_LINE_SIZE       2
#define VTK_TRIANGLE_SIZE   3
#define VTK_TETRA_SIZE      4

// Types of output files
#define GMSH_STYLE  1
#define FLOW_DATA_FILE 2
#define BOTH_OUTPUT 3

// Types of data, that could be written to output file
#define OUT_VECTOR_INT_SCA      1
#define OUT_VECTOR_INT_VEC      2
#define OUT_VECTOR_FLOAT_SCA    3
#define OUT_VECTOR_FLOAT_VEC    4
#define OUT_VECTOR_DOUBLE_SCA   5
#define OUT_VECTOR_DOUBLE_VEC   6
#define OUT_ARRAY_INT_SCA       7
#define OUT_ARRAY_FLOAT_SCA     8
#define OUT_ARRAY_DOUBLE_SCA    9

// TODO : i definice typu by mely byt uvnitr trid ktere tyto definice potrebuji
// opet ? verejne/nevrejne ...

/* Temporary structure for storing data */
typedef std::vector<double> ScalarFloatVector;
typedef std::vector< vector<double> > VectorFloatVector;

/* Temporary structure for storing data */
typedef struct OutScalar {
    ScalarFloatVector   *scalars;
    string              name;
    string              unit;
} OutScalars;

/* Temporary structure for storing data */
typedef struct OutVector {
    VectorFloatVector   *vectors;
    string              name;
    string              unit;
} OutVector;

/* Temporary structure for storing data */
typedef std::vector<OutScalar> OutScalarsVector;
typedef std::vector<OutVector> OutVectorsVector;

/**
 * Class of output data storing reference on data
 *
 * TODO: podrobnejsi komentar, metody komentovat v header filech
 * stejne tak u dlasich trid.
 */
class OutputData {
private:
public:
    // TODO: promenne by mely byt private

    string          *name;      ///< String with name of data
    string          *units;     ///< String with units
    void            *data;      ///< Pointer at own data
    unsigned char   type;       ///< Type values in vector
    int             comp_num;   ///< Number of components in vector
    int             num;        ///< Number of values in vector/array
    OutputData() {};            ///< Un-named constructor can't be called
    string* getName(void) { return name; };
    string* getUnits(void) { return units; };
    int getCompNum(void) { return comp_num; };
    int getValueNum(void) { return num; };
    OutputData(std::string name, std::string unit, int *data, unsigned int size);
    OutputData(std::string name, std::string unit, float *data, unsigned int size);
    OutputData(std::string name, std::string unit, double *data, unsigned int size);
    OutputData(std::string name, std::string unit, std::vector<int> &data);
    OutputData(std::string name, std::string unit, std::vector< vector<int> > &data);
    OutputData(std::string name, std::string unit, std::vector<float> &data);
    OutputData(std::string name, std::string unit, std::vector< vector<float> > &data);
    OutputData(std::string name, std::string unit, std::vector<double> &data);
    OutputData(std::string name, std::string unit, std::vector< vector<double> > &data);
    ~OutputData();
};

typedef std::vector<OutputData> OutputDataVec;

/**
 * Class of output
 */
class Output {
private:
    struct OutScalar *node_scalar;      // Temporary solution
    struct OutScalar *element_scalar;   // Temporary solution
    struct OutVector *element_vector;   // Temporary solution

    ofstream    *base_file;             ///< Base output stream
    string      *base_filename;         ///< Name of base output file
    string      *data_filename;         ///< Name of data output file
    ofstream    *data_file;             ///< Data output stream (could be same as base_file)
    int         format_type;            ///< Type of output
    Mesh        *mesh;
    std::vector<OutputData> *node_data; ///< List of data on nodes
    std::vector<OutputData> *elem_data; ///< List of data on elements


    // TODO: tohle by se melo resit jinak. napr: Output by melo obsahovat
    // abstraktni virtualni metodu write_data
    // pak bychom meli potomky tridy Output pro jednotlive formaty vystupu,
    // a ty by definovali vystup do patricneho formatu

    // Internal API for file formats
    int (*_write_data)(Output *output);
protected:
    // Protected getters for descendant

    // Protected setters for descendant
    void set_mesh(Mesh *_mesh) { mesh = _mesh; };
    void set_base_file(ofstream *_base_file) { base_file = _base_file; };
    void set_base_filename(string *_base_filename) { base_filename = _base_filename; };
    void set_format_type(int _format_type) { format_type = _format_type; };
    void set_node_data(std::vector<OutputData> *_node_data) { node_data = _node_data; };
    void set_elem_data(std::vector<OutputData> *_elem_data) { elem_data = _elem_data; };
public:

    // TODO: pokud je zde default konstruktor,
    // mela by byt moznost pozdeji nastavit ukazatel na sit, jinak se prislusna instance neda k nicemu pouzit

    Output() { node_scalar = NULL; element_scalar = NULL; element_vector = NULL; };

    Output(Mesh *mesh, string filename);
    ~Output();

    // Temporary solution
    void get_data_from_mesh(void);
    void free_data_from_mesh(void);

    template <typename _Data>
    int register_node_data(std::string name, std::string unit, _Data *data, uint size);

    template <typename _Data>
    int register_elem_data(std::string name, std::string unit, _Data *data, uint size);

    template <typename _Data>
    int register_node_data(std::string name, std::string unit, std::vector<_Data> &data);

    template <typename _Data>
    int register_elem_data(std::string name, std::string unit, std::vector<_Data> &data);

    // Public getters
    std::vector<OutputData> *get_node_data(void) { return node_data; };
    std::vector<OutputData> *get_elem_data(void) { return elem_data; };
    ofstream& get_base_file(void) { return *base_file; };
    string& get_base_filename(void) { return *base_filename; };
    ofstream& get_data_file(void) { return *data_file; };
    string& get_data_filename(void) { return *data_filename; };
    Mesh *get_mesh(void) { return mesh; };
    char get_format_type(void) { return format_type; };

    // Public setters
    void set_data_file(ofstream *_data_file) { data_file = _data_file; };

    int write_data(void);
};

/**
 * Class of output of during time
 */
class OutputTime : public Output {
private:
    int              current_step;      ///< Current step
    struct OutScalar *element_scalar;   // Temporary solution
    int              elem_sca_count;    // Temporary solution
    struct OutVector *element_vector;   // Temporary solution

    // Internal API for file formats
    int (*_write_data)(OutputTime *output, double time, int step);
    int (*_write_head)(OutputTime *output);
    int (*_write_tail)(OutputTime *output);
public:
    // Constructor and destructor
    OutputTime(Mesh *mesh, string filename);
    ~OutputTime();

    // Temporary solution for getting data from transport
    void get_data_from_transport(ConvectionTransport *transport);
    void free_data_from_transport(void);

    template <typename _Data>
    int register_node_data(std::string name, std::string unit, _Data *data, uint size);

    template <typename _Data>
    int register_elem_data(std::string name, std::string unit, _Data *data, uint size);

    // This method registers node data, that will be written to the file,
    // when write_data() will be called
    template <typename _Data>
    int register_node_data(std::string name, std::string unit, std::vector<_Data> &data);
    // This method register element data
    template <typename _Data>
    int register_elem_data(std::string name, std::string unit, std::vector<_Data> &data);

    // This method write data to the file
    int write_data(double time);
};

/* TODO: move to other file */
void output_flow_field_init(char *fname);
void output_flow_field_in_time(struct Problem *problem, double time);

/* TODO: move to new output_vtk.hh */
// Static data
int write_vtk_data(Output *output);
// Dynamic data
int write_vtk_head(OutputTime *output);
int write_vtk_time_data(OutputTime *output, double time, int step);
int write_vtk_tail(OutputTime *output);

/* TODO: move to new output_msh.hh */
// Static data
int write_msh_data(Output *output);
// Dynamic data
int write_msh_head(OutputTime *output);
int write_msh_time_data(OutputTime *output, double time, int step);
int write_msh_tail(OutputTime *output);

#endif
