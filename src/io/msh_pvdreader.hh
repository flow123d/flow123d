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
 * @file    msh_pvdreader.hh
 * @brief
 * @author  dalibor
 */

#ifndef MSH_PVD_READER_HH
#define	MSH_PVD_READER_HH

#include <string>
#include <istream>

#include "io/msh_basereader.hh"
#include "system/file_path.hh"


class VtkMeshReader;


class PvdMeshReader : public BaseMeshReader {
public:
	/**
     * Construct the PVD reader from given FilePath.
     * This opens the file for reading.
     */
	PvdMeshReader(const FilePath &file_name);

	/// Destructor
	~PvdMeshReader();

    /**
     * Read regions from the VTK file and save the physical sections as regions in the RegionDB.
     *
     * Region Labels starting with '!' are treated as boundary regions. Elements of these regions are used just to
     * assign regions to the boundary and are not used in actual FEM computations.
     */
    void read_physical_names(Mesh * mesh) override;

    /**
	 * Find header of DataArray section of VTK file given by field_name.
	 */
	MeshDataHeader & find_header(HeaderQuery &header_query) override;

protected:
	/// Represents data of one VTK file defined in PVD file.
	struct VtkFileData {
		/// Constructor
		VtkFileData(double t, FilePath file) : time(t), file_name(file), reader(nullptr) {}

		double time;             ///< time step of VTK file
		FilePath file_name;      ///< path to VTK file
		VtkMeshReader * reader;  ///< reader is created only if it's needed
	};

    /**
     * private method for reading of nodes
     */
    void read_nodes(Mesh * mesh) override;

    /**
     * Method for reading of elements.
     * Input of the mesh allows changing regions within the input file.
     */
    void read_elements(Mesh * mesh) override;

    /**
     * This method is specified for PVD reader. Table of mesh data headers (same as for GMSH or VTK) is not created,
     * but list of times and appropriate VTK files is filled.
     */
    void make_header_table() override;

    /**
     * Implements @p BaseMeshReader::read_element_data.
     */
    void read_element_data(ElementDataCacheBase &data_cache, MeshDataHeader header,
    		bool boundary_domain) override;

    /// Store list of VTK files and time steps declared in PVD file.
    std::vector<VtkFileData> file_list_;

    /// Path to PVD file allows construct FilePath objects of VTK files.
    std::string pvd_path_dir_;

    /// Iterator to items of \p file_list_
    std::vector<VtkFileData>::iterator list_it_;

};

#endif	/* MSH_PVD_READER_HH */

