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
#include "mesh/region.hh"
#include "mesh/elements.h"
#include "mesh/mesh.h"
#include "mesh/sides.h"
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
    ElementAccessor(const Mesh *mesh, unsigned int idx, bool boundary = false)
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
        return &(mesh_->element_vec_[element_idx_]);
    }
    
    inline arma::vec::fixed<spacedim> centre() const {
    	ASSERT(is_valid()).error("Invalid element accessor.");
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

    inline unsigned int index() const {
    	return (unsigned int)mesh_->elem_index(element_idx_);
    }

    inline void inc() {
        ASSERT(!is_regional()).error("Do not call inc() for regional accessor!");
        element_idx_++;
        r_idx_ = element()->region_idx();
        dim_=element()->dim();
    }

    inline SideIter side(const unsigned int loc_index) {
        return SideIter( Side(mesh_, element_idx_, loc_index) );
    }

    inline const SideIter side(const unsigned int loc_index) const {
        return SideIter( Side(mesh_, element_idx_, loc_index) );
    }

    bool operator==(const ElementAccessor<spacedim>& other) {
    	return (element_idx_ == other.element_idx_);
    }

    /**
     * -> dereference operator
     *
     * Allow simplify calling of element() method. Example:
 @code
     ElementAccessor<3> elm_ac(mesh, index);
     arma::vec centre;
     centre = elm_ac.element()->centre();  // full format of access to element
     centre = elm_ac->centre();            // short format with dereference operator
 @endcode
     */
    inline const Element * operator ->() const {
    	return &(mesh_->element_vec_[element_idx_]);
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
    /// True if the element is boundary
    bool boundary_;

    /// Index into Mesh::element_vec_ array.
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
