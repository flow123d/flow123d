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
 * @file    bc_mesh.h
 * @brief
 */

#ifndef BC_MESH_H
#define BC_MESH_H

#define NOT_IMPLEMENTED { ASSERT(false); return 0; }


#include "mesh/mesh.h"

class Partitioning;
class Neighbour;


/**
 * \brief Class represents boundary part of mesh.
 *
 * Holds pointer to parent Mesh and overwrites methods, that allow access to boundary
 * elements and other functionality and structures necessary to work above boundary
 * part of mesh.
 */
class BCMesh : public MeshBase {
public:
	/**
	 * Constructor from parent (bulk) Mesh.
	 */
	BCMesh(Mesh* parent_mesh);

	/// Destructor
	~BCMesh() override;

    /// Returns range of boundary elements of parent mesh
    Range<ElementAccessor<3>> elements_range() const override;

    /// Overwrite Mesh::get_part()
    Partitioning *get_part() override;

    /// Overwrite Mesh::get_local_part()
    const LongIdx *get_local_part() override;

    /// Overwrite Mesh::check_compatible_mesh()
    std::shared_ptr<std::vector<LongIdx>> check_compatible_mesh( Mesh & input_mesh) override;

    /// Overwrite Mesh::check_compatible_discont_mesh()
    std::shared_ptr<std::vector<LongIdx>> check_compatible_discont_mesh( Mesh & input_mesh) override;

    /// Overwrite Mesh::element_accessor()
    ElementAccessor<3> element_accessor(unsigned int idx) const override;

    /// Implement MeshBase::bc_mesh()
    BCMesh *bc_mesh() const override {
        return nullptr;
    }

private:

    /// setup distribution of elements and related vectors
    void init_distribution();

    // unused methods (should not be used)
    NodeAccessor<3> node(unsigned int) const override;
    Boundary boundary(unsigned int) const override;
    void check_element_size(unsigned int) const override;
    const RegionDB &region_db() const override;
    const DuplicateNodes *duplicate_nodes() const override NOT_IMPLEMENTED;



    /// Pointer to parent (bulk) mesh
	Mesh *parent_mesh_;

	/// Distribution of boundary elements to processors
	LongIdx *local_part_;


    friend class Mesh;

};


#endif  //BC_MESH_H
//-----------------------------------------------------------------------------
// vim: set cindent:
