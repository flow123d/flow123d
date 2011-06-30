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
 * @file    output_msh.h
 * @brief   Header: The functions for MSH (GMSH) outputs.
 *
 */

#ifndef OUTPUT_MSH_HH_
#define OUTPUT_MSH_HH_

/**
 * \brief This class is used for output data to VTK file format
 */
class OutputMSH : protected Output {
public:
    /**
     * \brief The constructor of this class
     */
    OutputMSH();

    /**
     * \brief The destructor of this class
     */
    ~OutputMSH();
private:
};

// TODO: make methods of OutputMSH from following functions

// Static data
int write_msh_data(Output *output);
// Dynamic data
int write_msh_head(OutputTime *output);
int write_msh_time_data(OutputTime *output, double time, int step);
int write_msh_tail(OutputTime *output);

#endif /* OUTPUT_MSH_HH_ */
