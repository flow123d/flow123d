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

#include "system.hh"

/// external types
struct Problem;
struct Element;
struct Transport;
class Mesh;

//=============================================================================
// TEMPORARY STRUCTURES
//=============================================================================
struct TTNode{
    double **conc;
    double *scalar;
};
struct TElement{
    double **conc;
    double *scalar;
    double **vector;
};
struct tripple {
    float d[3];
};
typedef std::vector<float> ScalarFloatVector;
typedef std::vector<tripple> VectorFloatVector;

/* TODO: convert to class */
typedef struct OutScalar {
    ScalarFloatVector *scalars;
    char name[32];
} OutScalars;

/* TODO: convert to class */
typedef struct OutVector {
    VectorFloatVector *vectors;
    char name[32];
} OutVector;

typedef std::vector<OutScalar> OutScalarsVector;
typedef std::vector<OutVector> OutVectorsVector;

// FILE formats
#define POS_ASCII           1
#define POS_BIN             2
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

// types of output files
#define GMSH_STYLE  1
#define FLOW_DATA_FILE 2
#define BOTH_OUTPUT 3

void output( struct Problem *problem );
void output_flow_field_init(struct Problem *problem);
void output_flow_field_in_time(struct Problem *problem,double time);
void output_init(struct Problem *problem);
void output_time(struct Problem *problem, double time);
FILE **open_temp_files(struct Transport *transport,const char *fileext,const char *open_param);

void output_msh_init_bin(Mesh*, char*);
void output_msh_init_ascii(Mesh*, char*);

void output_msh_init_vtk_serial_ascii( char *file);
void output_msh_finish_vtk_serial_ascii( char *file);
void output_transport_time_bin(struct Transport *transport, double time,int step,char *file);
void output_transport_time_ascii(struct Transport *transport, double time,int step,char *file);
void output_transport_time_vtk_serial_ascii(struct Transport *transport, double time, int step, char *file);
void write_ascii_header(struct Problem *problem, FILE *out);
void write_transport_ascii_data(FILE *out,struct Problem *problem,struct TTNode **nodes,struct TElement **elements,int time_steps,int ph);
void write_transport_binary_data(FILE *out,struct Problem *problem,struct TTNode **nodes,struct TElement **elements,int time_steps,int ph);

void output_flow_time_vtk_serial_ascii(Mesh *mesh,
        double time,
        int step,
        char *file);

#endif
//-----------------------------------------------------------------------------
// vim: set cindent:
