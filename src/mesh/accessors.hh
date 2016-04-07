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
 * @file    accessors.hh
 * @brief   
 */

#ifndef ACCESSORS_HH_
#define ACCESSORS_HH_

#include "mesh/bounding_box.hh"
#include "mesh/mesh_types.hh"
#include "mesh/region.hh"
#include "mesh/elements.h"
#include "mesh/mesh.h"
#include <armadillo>

/**
 * Element accessor templated just by dimension of the embedding space, used by Fields.
 * This should allow algorithms over elements where dimension of particular element is runtime parameter.
 *
 * This class suites as interface of Fields to the mesh elements, in particular this accessor knows directly
 * the region, and also can be used as an accessor that works on the whole region if used by Fields that do not depend on
 * particular elements as FieldConstant, FiledFormula, and FieldPython.
 *
 * TODO:
 * - make this kind of accessor subclass of FieldCommonBase or at least move it into src/fields
 *   since it has functionality particular for Fields
 *
 * Ideas:
 * need function to calculate intersection (object) of two ElementAccessors, but this definitely should be templated by
 * dimension of the ref. element (or rather shape of ref. element), here we can have case dispatch
 *
 */
template <int spacedim>
class ElementAccessor {
public:
    /**
     * Default invalid accessor.
     */
    ElementAccessor()
    : mesh_(NULL)
    {}

    /**
     * Regional accessor.
     */
    ElementAccessor(const Mesh *mesh, RegionIdx r_idx)
    : dim_(undefined_dim_), mesh_(mesh), r_idx_(r_idx)
    {}

    /**
     * Element accessor.
     */
    ElementAccessor(const Mesh *mesh, unsigned int idx, bool boundary)
    : mesh_(mesh), boundary_(boundary), element_idx_(idx), r_idx_(element()->region_idx())
    {
       dim_=element()->dim();
    }

    inline bool is_regional() const {
        return dim_ == undefined_dim_;
    }

    inline bool is_elemental() const {
        return ( is_valid() && ! is_regional() );
    }

    inline bool is_valid() const {
        return mesh_ != NULL;
    }

    inline unsigned int dim() const
        { return dim_; }

    inline const Element * element() const {
        if (boundary_) return (Element *)(mesh_->bc_elements(element_idx_)) ;
        else return  (Element *)(mesh_->element(element_idx_)) ;
    }

    inline arma::vec::fixed<spacedim> centre() const {
    	OLD_ASSERT(is_valid(), "Invalid element accessor.");
        if (is_regional() ) return arma::vec::fixed<spacedim>();
        else return element()->centre();
    }


    inline Region region() const
        { return Region( r_idx_, mesh_->region_db()); }

    inline RegionIdx region_idx() const
        { return r_idx_; }

    /// We need this method after replacing Region by RegionIdx, and movinf RegionDB instance into particular mesh
    //inline unsigned int region_id() const {
    //    return region().id();
    //}

    inline bool is_boundary() const {
        return boundary_;
    }

    inline unsigned int idx() const {
        return element_idx_;
    }
private:
    /**
     * When dim_ == undefined_dim_ ; the value of element_idx_ is invalid.
     * Is used for ElementAccessors for whole region
     */
    static const unsigned int undefined_dim_ = 100;

    /// Dimension of reference element.
    unsigned int dim_;

    /// Pointer to the mesh owning the element.
    const Mesh *mesh_;
    /// True if the element is boundary, i.e. stored in Mesh::bc_elements, bulk elements are stored in Mesh::element
    bool boundary_;

    /// Index into Mesh::bc_elements or Mesh::element array.
    unsigned int element_idx_;

    /// Region index.
    RegionIdx r_idx_;
};




/******************************************************************* implementations
 *
 *
 */
/*
template<int spacedim>
const BoundingBox &ElementAccessor<spacedim>::bounding_box() {
    return box_;
}
*/

#endif /* ACCESSORS_HH_ */
