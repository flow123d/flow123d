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
 * @file    output_vtk.h
 * @brief   Header: The functions for VTK outputs.
 */

#ifndef OUTPUT_VTK_HH_
#define OUTPUT_VTK_HH_

#include "input/accessors_forward.hh"

#include "output_data_base.hh"
#include "output_time.hh"

#include <ostream>

using namespace std;

/**
 * \brief This class is used for output data to VTK file format
 */
class OutputVTK : public OutputTime {

public:
	typedef OutputTime FactoryBaseType;

    /**
     * \brief The constructor of this class. The head of file is written, when
     * constructor is called
     */
    OutputVTK(const Input::Record &in_rec);

    /**
     * \brief The constructor of this class. The head of file is written, when
     * constructor is called
     */
    OutputVTK();

    /**
     * \brief The destructor of this class. It writes tail of the file too.
     */
    ~OutputVTK();

    /**
     * \brief The definition of input record for vtk file format
     */
    static const Input::Type::Record & get_input_type();

    /**
	 * \brief The definition of input record for selection of variant of file format
	 */
    static const Input::Type::Selection & get_input_type_variant();

    /**
	 * \brief The definition of input record for selection of compression type
	 */
    static const Input::Type::Selection & get_input_type_compression();

    /**
     * \brief This function write data to VTK (.pvd) file format
     * for curent time
     */
    int write_data(void);

    /**
     * \brief This function writes header of VTK (.pvd) file format
     */
    int write_head(void);

    /**
     * \brief This function writes tail of VTK (.pvd) file format
     */
    int write_tail(void);

protected:

    /**
     * \brief The declaration enumeration used for variant of file VTK format
     */
    typedef enum Variant {
    	VARIANT_ASCII  = 1,
    	VARIANT_BINARY = 2
    } Variant;

    /**
     * \brief The declaration of enumeration used for type of compression
     * used in file format
     */
    typedef enum Compression {
    	COMPRESSION_NONE = 1,
    	COMPRESSION_GZIP = 2
    } Compression;

    // VTK Element types
    typedef enum {
        VTK_VERTEX = 1,
        VTK_POLY_VERTEX = 2,
        VTK_LINE = 3,
        VTK_POLY_LINE = 4,
        VTK_TRIANGLE = 5,
        VTK_TRIANGLE_STRIP = 6,
        VTK_POLYGON = 7,
        VTK_PIXEL = 8,
        VTK_QUAD = 9,
        VTK_TETRA = 10,
        VTK_VOXEL = 11,
        VTK_HEXAHEDRON = 12,
        VTK_WEDGE = 13,
        VTK_PYRAMID = 14,
        VTK_QUADRIC_EDGE = 21,
        VTK_QUADRIC_TRIANGLE = 22,
        VTK_QUADRIC_QUAD = 23,
        VTK_QUADRIC_TETRA = 24,
        VTK_QUADRIC_HEXAHEDRON = 25
    } VTKElemType;

    // VTK Element size (number of nodes)
    typedef enum {
        VTK_LINE_SIZE = 2,
        VTK_TRIANGLE_SIZE = 3,
        VTK_TETRA_SIZE = 4
    } VTKElemSize;

    typedef enum { VTK_INT8, VTK_UINT8, VTK_INT16, VTK_UINT16, VTK_INT32, VTK_UINT32, 
                   VTK_FLOAT32, VTK_FLOAT64
    } VTKValueType;

    static const std::string vtk_value_type_map(VTKValueType t) {
        static const std::vector<std::string> types = {
            "Int8", "UInt8", "Int16", "UInt16", "Int32", "UInt32",
            "Float32","Float64"};
        return types[t];
    };

    /// Registrar of class to factory
    static const int registrar;

    /**
     * \brief Write header of VTK file (.vtu)
     */
    void write_vtk_vtu_head(void);

    /**
     * \brief Fills the given vector with VTK element types indicators.
     */
    void fill_element_types_vector(std::vector<unsigned int> &data);

    /**
     * Write registered data to output stream using ascii format
     */
    void write_vtk_data_ascii(OutputDataFieldVec &output_data_map);

    /**
     * Write registered data to output stream using ascii format
     */
    void write_vtk_data_ascii(OutputDataPtr output_data, VTKValueType type);
    
    /**
     * \brief Write names of data sets in @p output_data vector that have value type equal to @p type.
     * Output is done into stream @p file.
     */
    void write_vtk_data_names(ofstream &file,
            OutputDataFieldVec &output_data_map);

    /**
     * \brief Write data on nodes to the VTK file (.vtu)
     */
    void write_vtk_node_data(void);

    /**
     * \brief Write data on elements to the VTK file (.vtu)
     */
   void write_vtk_element_data(void);

   /**
    * \brief Write tail of VTK file (.vtu)
    */
   void write_vtk_vtu_tail(void);

   /**
    * \brief This function write all scalar and vector data on nodes and elements
    * to the VTK file (.vtu)
    */
   void write_vtk_vtu(void);

   /**
    * Set appropriate file path substrings.
    * Make subdirectory for VTU time frames.
    */
   void make_subdirectory();


   /**
    * Data output stream (could be same as base_file)
    */
   ofstream _data_file;

   /**
    * Path to time frame VTU data subdirectory
    */
   string subdir_name_;

   string main_output_basename_;

   string main_output_dir_;
};

#endif /* OUTPUT_VTK_HH_ */
