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

#include "io/output.h"

/**
 * \brief This class is used for output data to VTK file format
 */
class OutputMSH : public OutputTime {
public:

    /**
     * \brief The constructor of this class
     */
    OutputMSH();

    /**
     * \brief The constructor of this class
     */
    OutputMSH(const Input::Record &in_rec);

    /**
     * \brief The destructor of this class
     */
    ~OutputMSH();

    /**
     * \brief The definition of input record for gmsh file format
     */
    static Input::Type::Record input_type;

    /**
     * \brief This method writes head of GMSH (.msh) file format
     *
     * \return      This function returns 1
     */
    int write_head(void);

    /**
     * \brief This method writes data to GMSH (.msh) file format for current time
     *
     * \return      This function returns 1
     */
    int write_data(void);

    /**
     * \brief This method should write tail of GMSH (.msh) file format
     *
     * It is stupid file format. It doesn't write anything special at the end of
     * the file
     *
     * \return      This function returns 1
     */
    int write_tail(void);


private:

    /**
     * Was header already written to output file?
     */
    bool header_written;

    /**
     * \brief This function write header of GMSH (.msh) file format
     */
    void write_msh_header(void);

    /**
     * \brief This function writes geometry (position of nodes) to GMSH (.msh) file
     * format
     */
    void write_msh_geometry(void);

    /**
     * \brief This function writes topology (connection of nodes) to the GMSH (.msh)
     * file format
     */
    void write_msh_topology(void);

    /**
     * \brief This function writes continuous ascii data to GMSH (.msh) output file.
     *
     * \param[in]   *out_data   The pointer at structure storing pointer at own data.
     */
    void write_msh_ascii_cont_data(OutputDataBase* output_data);

    /**
     * \brief This function writes discontinuous ascii data to GMSH (.msh) output file.
     *
     * \param[in]   *out_data   The pointer at structure storing pointer at own data.
     */
    void write_msh_ascii_discont_data(OutputDataBase* output_data);

    /**
     * \brief This function write all data on nodes to output file. This function
     * is used for static and dynamic data
     *
     * \param[in]   time        The time from start
     * \param[in]   step        The number of steps from start
     */
    void write_msh_node_data(double time, int step);


    /**
     * \brief This function write all data on elements to output file. This
     * function is used for static and dynamic data
     *
     * \param[in]   time        The time from start
     * \param[in]   step        The number of steps from start
     */
    void write_msh_elem_data(double time, int step);
};

#endif /* OUTPUT_MSH_HH_ */
