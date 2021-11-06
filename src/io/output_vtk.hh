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

#include <memory>          // for shared_ptr
#include <ostream>         // for ofstream, stringstream, ostringstream
#include <string>          // for string
#include "output_time.hh"  // for OutputTime, OutputTime::OutputDataFieldVec
#include <zlib.h>

class Mesh;
class TimeUnitConversion;
namespace Input {
	class Record;
	namespace Type {
		class Record;
		class Selection;
	}
}
template <typename T> class ElementDataCache;

using namespace std;


/**
 * \brief This class is used for output data to VTK file format
 */
class OutputVTK : public OutputTime {

public:
	typedef OutputTime FactoryBaseType;
	typedef uLong zlib_ulong;


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

    /// Override @p OutputTime::init_from_input.
    void init_from_input(const std::string &equation_name,
                         const Input::Record &in_rec,
                         const std::shared_ptr<TimeUnitConversion>& time_unit_conv) override;

protected:

    /**
     * \brief The declaration enumeration used for variant of file VTK format
     */
    typedef enum {
    	VARIANT_ASCII  = 0,
    	VARIANT_BINARY_UNCOMPRESSED = 1,
    	VARIANT_BINARY_ZLIB = 2
    } VTKVariant;

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

    /// Registrar of class to factory
    static const int registrar;

    /// Formats of DataArray section
	static const std::vector<std::string> formats;

	/**
	 * Used internally by write_data.
	 */
	string form_vtu_filename_(string basename, int i_step, int rank);

	/**
     * \brief Write header of VTK file (.vtu)
     */
    void write_vtk_vtu_head(void);

    /**
     * \brief Fills the data cache with VTK element types indicators.
     */
    std::shared_ptr<ElementDataCache<unsigned int>> fill_element_types_data();

    /**
     * Write registered data of all components of given Field to output stream
     */
    void write_vtk_field_data(OutputDataFieldVec &output_data_map);

    /**
     * Write output data stored in OutputData vector to output stream
     */
    void write_vtk_data(OutputDataPtr output_data, unsigned int start = 0);
    
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
    * \brief Write native data (part of our own data skipped by Paraview) to the VTK file (.vtu)
    *
    * Tags of native data are subtags of 'Flow123dData' tag, that is subtag of 'Piece' tag
    */
  void write_vtk_native_data(void);

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
    * Compress data stored in @p uncompressed_stream to @p compressed_stream.
    *
    * Use ZLib compression.
    */
   void compress_data(stringstream &uncompressed_stream, stringstream &compressed_stream);


   /**
    * Data output stream (could be same as base_file)
    */
   ofstream _data_file;

   /**
    * Stream of appended data (used only for binary appended output)
    */
   ostringstream appended_data_;

   /**
    * Path to time frame VTU data subdirectory
    */
   string subdir_name_;

   /// Basename of main output file (without extension)
   string main_output_basename_;

   /// Main output file directory
   string main_output_dir_;

   /// Output format (ascii, binary or binary compressed)
   VTKVariant variant_type_;
};

#endif /* OUTPUT_VTK_HH_ */
