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
 * @file    raw_mesh.hh
 * @ingroup mesh
 * @brief   Mesh construction
 */

#ifndef RAW_MESH_HH_
#define RAW_MESH_HH_

#include <mpi.h>                             // for MPI_Comm, MPI_COMM_WORLD
#include "mesh/elements.h"
#include "input/accessors.hh"                // for Record, Array (ptr only)
#include "system/armor.hh"

class Mesh;
class RegionIdx;
class Distribution;

class RawMesh {
public:
    /**
     * Empty constructor.
     *
     * Use only for unit tests!!!
     */
	RawMesh();
    /**
     * Constructor from mesh.
     * Do not process input record. That is done in init_from_input.
     */
	RawMesh(Mesh *mesh);

	/// Copy constructor
	RawMesh(RawMesh &other);

    /// Destructor.
    ~RawMesh();

private:
    /**
     * Part of the constructor whichdoes not depedn on input record.
     * Initializes node-side numbering according to RefElement.
     */
    void init(Mesh *mesh);

    void add_element(Mesh *mesh, unsigned int dim, RegionIdx region_idx, unsigned int partition_id, std::vector<unsigned int> node_ids);

    /**
     * Accessor to the input record for the mesh.
     */
    Input::Record in_record_;

    /**
     * MPI communicator used for partitioning and ...
     */
    MPI_Comm comm_;

    /**
     * Vector of nodes of the mesh.
     */
    shared_ptr<Armor::Array<double>> nodes_;

    /**
     * Vector of elements of the mesh.
     */
    vector<Element> element_vec_;

    /// Partition numbers for local elements in original distribution of elements given be @p init_el_ds_.
    LongIdx * loc_part_;

    /// Original distribution of elements. Depends on type of partitioner
    Distribution * init_el_ds_;
};


#endif // RAW_MESH_HH_
