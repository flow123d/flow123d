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

#include "io/output.h"

/**
 * \brief This class is used for output data to VTK file format
 */
class OutputVTK : public OutputFormat {

public:
    /**
     * \brief The constructor of this class. The head of file is written, when
     * constructor is called
     */
    OutputVTK(Output *_output);

    /**
     * \brief The constructor of this class. The head of file is written, when
     * constructor is called
     */
    OutputVTK(OutputTime *_output_time, const Input::Record &in_rec);

    /**
     * \brief The constructor of this class. The head of file is written, when
     * constructor is called
     */
    OutputVTK(OutputTime *_output_time);

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
     * \brief This function output data to serial VTK file format (single .vtu)
     */
    int write_data(void);

    /**
     * \brief This function write data to VTK (.pvd and .vtu) file format
     * for specific time
     *
     * \param[in]   time        The time from start
     */
    int write_data(double time);

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

    /**
     * \brief The pointer at Output
     */
    Output *output;

    /**
     * \brief The pointer at OutputTime
     */
    OutputTime *output_time;

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
     * \brief This function writes ascii continuous data to VTK (.vtu) output file.
     *
     * \param[in]   *data   The pointer at structure storing pointer at own data.
     */
    void write_vtk_ascii_cont_data(OutputData *data);

    /**
     * \brief This function writes ascii continuous data to discontinuous mesh
     * to VTK (.vtu) output file.
     *
     * \param[in]   *data   The pointer at structure storing pointer at own data.
     */
    void write_vtk_ascii_discont_data(OutputData *data);

    /**
     * \brief This function writes ascii data to VTK (.vtu) output file.
     *
     * \param[in]   *data   The pointer at structure storing pointer at own data.
     */
    void write_vtk_ascii_data(OutputData *data);

    /**
     * \brief Write scalar data to the VTK file (.vtu)
     *
     * \param[in]   *data   The pointer at structure storing pointer at own data.
     */
    void write_vtk_scalar_ascii(OutputData *data);

    /**
     * \brief Write vector data to VTK file (.vtu)
     *
     * \param[in]   *data   The pointer at structure storing pointer at own data.
     */
    void write_vtk_vector_ascii(OutputData *data);

    /**
     * \brief Go through all vectors of scalars and vectors and call functions that
     * write these data to VTK file (.vtu)
     *
     * \param[in]   *data   The pointer at vector of data
     */
    void write_vtk_data_ascii(std::vector<OutputData> *data);

    /**
     * \brief Write names of scalar values to the VTK file (.vtu)
     *
     * \param[in]   *data   The pointer at vector of data
     */
    void write_vtk_scalar_data_names(vector<OutputData> *data);

    /**
     * \brief Write names of vector values to the VTK file (.vtu)
     *
     * \param[in]   *data   The pointer at vector of data
     */
    void write_vtk_vector_data_names(vector<OutputData> *data);

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
