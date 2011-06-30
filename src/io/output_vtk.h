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

#include "io/output.h"

// TODO: move enums to VTK class

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
 * \brief This class is used for output data to VTK file format
 */
class OutputVTK : protected Output {
public:
    /**
     * \brief The constructor of this class
     */
    OutputVTK();

    /**
     * \brief The destructor of this class
     */
    ~OutputVTK();
private:
};

// TODO: make methods of OutputVTK from following functions

// Static data
int write_vtk_vtu_data(Output *output);

// Dynamic data
int write_vtk_pvd_head(OutputTime *output);
int write_vtk_pvd_data(OutputTime *output, double time, int step);
int write_vtk_pvd_tail(OutputTime *output);

#endif /* OUTPUT_VTK_HH_ */
