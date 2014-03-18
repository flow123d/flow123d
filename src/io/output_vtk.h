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
 * $Id: $
 * $Revision: $
 * $LastChangedBy: $
 * $LastChangedDate: $
 *
 * @file    output_vtk.h
 * @brief   Header: The functions for VTK outputs.
 *
 */

#ifndef OUTPUT_VTK_HH_
#define OUTPUT_VTK_HH_

#include "input/accessors.hh"
#include "fields/field_base.hh"

#include "io/output.h"

/**
 * \brief This class is used for output data to VTK file format
 */
class OutputVTK : public OutputTime {

public:

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
    static Input::Type::Record input_type;

    /**
	 * \brief The definition of input record for selection of variant of file format
	 */
    static Input::Type::Selection input_type_variant;

    /**
	 * \brief The definition of input record for selection of compression type
	 */
    static Input::Type::Selection input_type_compression;

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

private:

    /**
     * Was header already written to output file?
     */
    bool header_written;

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

    /**
     * \brief Write header of VTK file (.vtu)
     */
    void write_vtk_vtu_head(void);

    /**
     * \brief Write geometry (position of nodes) to the VTK file (.vtu)
     */
    void write_vtk_geometry(void);

    /**
     * \brief Write topology (connection of nodes) to the VTK file (.vtu)
     */
    void write_vtk_topology(void);

    /**
     * \brief Write geometry (position of nodes) to the VTK file (.vtu)
     *
     * This method is used, when discontinuous data are saved to the .vtu file
     */
    void write_vtk_discont_geometry(void);

    /**
     * \brief Write topology (connection of nodes) to the VTK file (.vtu)
     *
     * This method is used, when discontinuous data are saved to the .vtu file
     */
    void write_vtk_discont_topology(void);

    /**
     *
     */
    void write_vtk_ascii_data(OutputDataBase *output_data);

    /**
     *
     */
    void write_vtk_data_ascii(vector<OutputDataBase*> &output_data);

    /**
     * \brief Write names of data sets in @p output_data vector that have value type equal to @p type.
     * Output is done into stream @p file.
     */
    void write_vtk_data_names(ofstream &file, vector<OutputDataBase*> &output_data);

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
};

#endif /* OUTPUT_VTK_HH_ */
