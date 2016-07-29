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
 * @file    output_msh.h
 * @brief   Header: The functions for MSH (GMSH) outputs.
 */

#ifndef OUTPUT_MSH_HH_
#define OUTPUT_MSH_HH_

#include "output_time.hh"
#include "system/sys_vector.hh"

/**
 * \brief This class is used for output data to VTK file format
 */
class OutputMSH : public OutputTime {
public:
	typedef OutputTime FactoryBaseType;

    /**
     * \brief The constructor of this class.
     * We open the output file in first call of write_data
     */
    OutputMSH();

    /**
     * \brief The destructor of this class
     */
    ~OutputMSH();

    /**
     * \brief The definition of input record for gmsh file format
     */
    static const Input::Type::Record & get_input_type();

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

    /// Registrar of class to factory
    static const int registrar;

    bool header_written;

    /**
     * EquationOutput force output of all output fields at the time zero.
     * We keep this list to perform single elemnent/node output in order to have correct behavior in GMSH.
     */
    std::vector< std::vector< OutputDataPtr >> dummy_data_list_;

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
    template<class element>
    void write_msh_ascii_cont_data(flow::VectorId<element> &vec, OutputDataPtr output_data);

    /**
     * \brief This function writes discontinuous ascii data to GMSH (.msh) output file.
     *
     * \param[in]   *out_data   The pointer at structure storing pointer at own data.
     */
    void write_msh_ascii_discont_data(OutputDataPtr output_data);

    /**
     * \brief This function write all data on nodes to output file. This function
     * is used for static and dynamic data
     *
     * \param[in]   time        The time from start
     * \param[in]   step        The number of steps from start
     */
    void write_node_data(OutputDataPtr output_data);
    /**
     * \brief writes ElementNode data ascii GMSH (.msh) output file.
     *
     */
    void write_corner_data(OutputDataPtr output_data);


    /**
     * \brief This function write all data on elements to output file. This
     * function is used for static and dynamic data
     *
     * \param[in]   time        The time from start
     * \param[in]   step        The number of steps from start
     */
    void write_elem_data(OutputDataPtr output_data);

    /**
     * Write fields of single discrete space type, use given format function for the output of single field data.
     *
     * At first call it fills dummy_data_list_ (assuming output of all fields at the first output frame).
     */
    void write_field_data(OutputTime::DiscreteSpace type_idx, void (OutputMSH::* format_fce)(OutputDataPtr) );

    /**
     * \brief This method add right suffix to .msh GMSH file
     */
    void fix_base_file_name(void);
};

#endif /* OUTPUT_MSH_HH_ */
