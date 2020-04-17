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

#include <vector>          // for vector
#include "output_time.hh"  // for OutputTime::OutputDataPtr, OutputTime, Out...
namespace Input { namespace Type { class Record; } }
namespace flow { template <class T> class VectorId; }


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

    /// Complete information about dummy fields that are not in output_data_list_.
    void add_dummy_fields() override;

    /**
     * Set shared pointers of output data caches.
     */
    void set_output_data_caches(std::shared_ptr<OutputMeshBase> mesh_ptr) override;

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
     * \brief This function writes ascii data to GMSH (.msh) output file.
     *
     * \param[in]   id_cache     Data cache of node or element ids.
     * \param[in]   output_data  The pointer at structure storing pointer at own data.
     * \param[in]   discont      Flag determines continuous or discontinuous mesh.
     */
    void write_msh_ascii_data(std::shared_ptr<ElementDataCache<unsigned int>> id_cache, OutputDataPtr output_data, bool discont = false);

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
     * \brief This method add right suffix to .msh GMSH file
     */
    void fix_base_file_name(void);

    /// Vector gets ids of nodes.
    std::shared_ptr<ElementDataCache<unsigned int>> node_ids_;
    /// Vector gets ids of elements.
    std::shared_ptr<ElementDataCache<unsigned int>> elem_ids_;
    /// Vector gets ids of regions.
    std::shared_ptr<ElementDataCache<unsigned int>> region_ids_;
    /// Vector gets partitions of elements.
    std::shared_ptr<ElementDataCache<int>> partitions_;
};

#endif /* OUTPUT_MSH_HH_ */
