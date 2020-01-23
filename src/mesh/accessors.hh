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
#include "mesh/node_accessor.hh"
#include "mesh/sides.h"
#include "la/distribution.hh"
#include "mesh/point.hh"
#include <vector>
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
    typedef typename Space<spacedim>::Point Point;
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
    ElementAccessor(const Mesh *mesh, unsigned int idx)
    : mesh_(mesh), boundary_(idx>=mesh->n_elements()), element_idx_(idx), r_idx_(element()->region_idx())
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

    /// Return local idx of element in boundary / bulk part of element vector
    inline unsigned int idx() const {
        if (boundary_) return ( element_idx_ - mesh_->bulk_size_ );
        else return element_idx_;
    }

    /// Return global idx of element in full element vector
    inline unsigned int mesh_idx() const {
        return element_idx_;
    }

    inline unsigned int index() const {
    	return (unsigned int)mesh_->find_elem_id(element_idx_);
    }
    
    inline unsigned int proc() const {
        return mesh_->get_el_ds()->get_proc(mesh_->get_row_4_el()[element_idx_]);
    }

    inline void inc() {
        ASSERT(!is_regional()).error("Do not call inc() for regional accessor!");
        element_idx_++;
        r_idx_ = element()->region_idx();
        dim_=element()->dim();
        boundary_ = (element_idx_>=mesh_->n_elements());
    }

    inline SideIter side(const unsigned int loc_index) {
        return SideIter( Side(mesh_, element_idx_, loc_index) );
    }

    inline const SideIter side(const unsigned int loc_index) const {
        return SideIter( Side(mesh_, element_idx_, loc_index) );
    }

//    inline const NodeAccessor<spacedim> node_accessor(unsigned int ni) {
//
//    }

    inline NodeAccessor<spacedim> node(unsigned int ni) const {
        return mesh_->node( element()->node_idx(ni) );
    }

    /**
    * Return bounding box of the element.
    * Simpler code, but need to check performance penelty.
    */
    inline BoundingBox bounding_box() const {
        return BoundingBox(this->vertex_list());
    }

    /**
     * Return list of element vertices.
     * TODO: remove, SLOW!!
     */
    inline vector<Point> vertex_list() const {
        vector<Point> vertices(element()->n_nodes());
        for(unsigned int i=0; i<element()->n_nodes(); i++)
            vertices[i] = *node(i);
        return vertices;
    }

    /// Computes the measure of the element.
    double measure() const;

    /** Computes the Jacobian of the element.
     * J = det ( 1  1  1  1 )
     *           x1 x2 x3 x4
     *           y1 y2 y3 y4
     *           z1 z2 z3 z4
     */
    inline double tetrahedron_jacobian() const
    {
        ASSERT(dim() == 3)(dim()).error("Cannot provide Jacobian for dimension other than 3.");
        return arma::dot( arma::cross(
                *node(1) - *node(0),
                *node(2) - *node(0)),
                *node(3) - *node(0)
                );
    }

    /// Computes the barycenter.
    Point centre() const;

    /**
     * Quality of the element based on the smooth and scale-invariant quality measures proposed in:
     * J. R. Schewchuk: What is a Good Linear Element?
     *
     * We scale the measure so that is gives value 1 for regular elements. Line 1d elements
     * have always quality 1.
     */
    double quality_measure_smooth(SideIter side) const;

    inline bool operator==(const ElementAccessor<spacedim>& other) const {
    	return (element_idx_ == other.element_idx_);
    }

    inline bool operator!=(const ElementAccessor<spacedim>& other) const {
    	return (element_idx_ != other.element_idx_);
    }

    /**
     * -> dereference operator
     *
     * Allow simplify calling of element() method. Example:
 @code
     ElementAccessor<3> elm_ac(mesh, index);
     arma::vec centre;
     centre = elm_ac.element()->node_idx(0);  // full format of access to element
     centre = elm_ac->node_idx(0);            // short format with dereference operator
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

/**
 * SET THE "METRICS" FIELD IN STRUCT ELEMENT
 */
template <int spacedim>
double ElementAccessor<spacedim>::measure() const {
    switch (dim()) {
        case 0:
            return 1.0;
            break;
        case 1:
            return arma::norm(*node(1) - *node(0), 2);
            break;
        case 2:
            return
                arma::norm(
                    arma::cross(*node(1) - *node(0), *node(2) - *node(0)),
                    2
                ) / 2.0 ;
            break;
        case 3:
            return fabs(
                arma::dot(
                    arma::cross(*node(1) - *node(0), *node(2) - *node(0)),
                    *node(3) - *node(0))
                ) / 6.0;
            break;
    }
    return 1.0;
}

/**
 * SET THE "CENTRE[]" FIELD IN STRUCT ELEMENT
 */

template <int spacedim>
typename Space<spacedim>::Point ElementAccessor<spacedim>::centre() const {
	ASSERT(is_valid()).error("Invalid element accessor.");
    if (is_regional() ) return Point();

    Point centre;
    centre.zeros();

    for (unsigned int li=0; li<element()->n_nodes(); li++) {
        centre += *node(li);
    }
    centre /= (double) element()->n_nodes();
    return centre;
}


template <int spacedim>
double ElementAccessor<spacedim>::quality_measure_smooth(SideIter side) const {
    if (dim_==3) {
        double sum_faces=0;
        double face[4];
        for(unsigned int i=0; i<4; i++, ++side) sum_faces+=( face[i]=side->measure());

        double sum_pairs=0;
        for(unsigned int i=0;i<3;i++)
            for(unsigned int j=i+1;j<4;j++) {
                unsigned int i_line = RefElement<3>::line_between_faces(i,j);
                auto line_node_1 = RefElement<3>::interact(Interaction<0,1>(i_line))[1];
                auto line_node_0 = RefElement<3>::interact(Interaction<0,1>(i_line))[0];
                arma::vec3 line = *node(line_node_1) - *node(line_node_0);
                sum_pairs += face[i]*face[j]*arma::dot(line, line);
            }
        double regular = (2.0*sqrt(2.0/3.0)/9.0); // regular tetrahedron
        return fabs( measure()
                * pow( sum_faces/sum_pairs, 3.0/4.0))/ regular;

    }
    if (dim_==2) {
        return fabs(
                measure()/
                pow(
                         arma::norm(*node(1) - *node(0), 2)
                        *arma::norm(*node(2) - *node(1), 2)
                        *arma::norm(*node(0) - *node(2), 2)
                        , 2.0/3.0)
               ) / ( sqrt(3.0) / 4.0 ); // regular triangle
    }
    return 1.0;
}

#endif /* ACCESSORS_HH_ */
