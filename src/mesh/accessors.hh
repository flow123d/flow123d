/*!
 *
 * Copyright (C) 2015 Technical University of Liberec.  All rights reserved.
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
#include "mesh/ref_element.hh"
#include "la/distribution.hh"
#include "mesh/point.hh"
#include <vector>
#include <armadillo>


/**
 * Due to circular dependence of return parameters in mesh accessors,
 * and the intention to have all methods inlined,
 * we gathered the accessors into a single file.
 * This way, it is possible to implement it, see the simple test below.
 * 
 * Do not include accessors_impl.hh file anywhere but at the end of this file!
 * 
 * previous include loops:
 * side -> boundary -> elemen accessor
 * element accessor <-> side
 * edge <-> side
 * boundary -> edge
 */

// Compilable test of class loop dependence of return parameters.
// class A;
// class B;

// class A{
//     public:
//       A()
//       {};
      
//       B create_b();
// };

// class B{
//     public:
//       B()
//       {};
      
//       A create_a();
// };


// B A::create_b()
// { return B(); }

// A B::create_a()
// { return A(); }


class Side;
class SiderIter;
class Edge;
class Boundary;

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
    /// Default invalid accessor.
    ElementAccessor();

    /// Regional accessor.
    ElementAccessor(const Mesh *mesh, RegionIdx r_idx);

    /// Element accessor.
    ElementAccessor(const Mesh *mesh, unsigned int idx);

    /// Incremental function of the Element iterator.
    void inc();

    /// Return list of element vertices.
    vector<arma::vec3> vertex_list() const;

    /// Computes the measure of the element.
    double measure() const;

    /** Computes the Jacobian of the element.
     * J = det ( 1  1  1  1 )
     *           x1 x2 x3 x4
     *           y1 y2 y3 y4
     *           z1 z2 z3 z4
     */
    double tetrahedron_jacobian() const;

    /// Computes the barycenter.
    arma::vec::fixed<spacedim> centre() const;

    /**
     * Quality of the element based on the smooth and scale-invariant quality measures proposed in:
     * J. R. Schewchuk: What is a Good Linear Element?
     *
     * We scale the measure so that is gives value 1 for regular elements. Line 1d elements
     * have always quality 1.
     */
    double quality_measure_smooth(SideIter side) const;

    SideIter side(const unsigned int loc_index);

    const SideIter side(const unsigned int loc_index) const;



    bool is_regional() const {
        return dim_ == undefined_dim_;
    }

    bool is_elemental() const {
        return ( is_valid() && ! is_regional() );
    }

    bool is_valid() const {
        return mesh_ != NULL;
    }

    unsigned int dim() const
        { return dim_; }

    const Element * element() const {
        return &(mesh_->element_vec_[element_idx_]);
    }
    

    Region region() const
        { return Region( r_idx_, mesh_->region_db()); }

    RegionIdx region_idx() const
        { return r_idx_; }

    /// We need this method after replacing Region by RegionIdx, and movinf RegionDB instance into particular mesh
    //unsigned int region_id() const {
    //    return region().id();
    //}

    bool is_boundary() const {
        return boundary_;
    }

    /// Return local idx of element in boundary / bulk part of element vector
    unsigned int idx() const {
        if (boundary_) return ( element_idx_ - mesh_->bulk_size_ );
        else return element_idx_;
    }

    /// Return global idx of element in full element vector
    unsigned int mesh_idx() const {
        return element_idx_;
    }

    unsigned int index() const {
    	return (unsigned int)mesh_->find_elem_id(element_idx_);
    }
    
    unsigned int proc() const {
        return mesh_->get_el_ds()->get_proc(mesh_->get_row_4_el()[element_idx_]);
    }


    NodeAccessor<3> node(unsigned int ni) const {
        return mesh_->node( element()->node_idx(ni) );
    }

    /**
    * Return bounding box of the element.
    * Simpler code, but need to check performance penelty.
    */
    BoundingBox bounding_box() const {
        return BoundingBox(this->vertex_list());
    }

    bool operator==(const ElementAccessor<spacedim>& other) const {
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
    const Element * operator ->() const {
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





//=============================================================================
// Edge class
//=============================================================================
class Edge
{
public:
    /// Default invalid edge accessor constructor.
    Edge();

    /// Valid edge accessor constructor.
    Edge(const Mesh *mesh, unsigned int edge_idx);

    /// Gets side iterator of the @p i -th side.
    SideIter side(const unsigned int i) const;



    bool is_valid() const
    { return mesh_ != nullptr; }

    /// Returns edge global index.
    unsigned int idx() const {
        ASSERT_DBG(is_valid());
        return edge_idx_;
    }

    /// Incremental function of the Edge iterator.
    void inc() {
        ASSERT_DBG(is_valid()).error("Do not call inc() for invalid accessor!");
        edge_idx_++;
    }

    /// Comparison operator of the iterator.
    bool operator==(const Edge& other) const{
    	return (edge_idx_ == other.edge_idx_);
    }

    /// Returns number of sides aligned with the edge.
    unsigned int n_sides() const
    { return edge_data()->n_sides;}

private:
    /// Pointer to the mesh owning the node.
    const Mesh *mesh_;
    /// Index into Mesh::edges vector.
    unsigned int edge_idx_;

    /// Getter for edge data from mesh.
    const EdgeData* edge_data() const;
};





//=============================================================================
// Boundary class
//=============================================================================
class Boundary
{
public:
    Boundary();
    Boundary(BoundaryData* boundary_data);

    Edge edge();
    ElementAccessor<3> element_accessor();
    Region region();
    Element * element();

    bool is_valid() const {
        return boundary_data_ != nullptr;
    }
    
    Mesh* mesh() {
        ASSERT_DBG(is_valid());
        return boundary_data_->mesh_;
    }

    uint edge_idx() {
        ASSERT_DBG(is_valid());
        return boundary_data_->edge_idx_;
    }

    uint bc_ele_idx() {
        ASSERT_DBG(is_valid());
        return boundary_data_->bc_ele_idx_;
    }

private:
    BoundaryData* boundary_data_;
};





//=============================================================================
// Side class
//=============================================================================
class Side {
public:
    /// Default invalid side accessor constructor.
    Side();

    /// Valid edge accessor constructor.
    Side(const Mesh * mesh, unsigned int elem_idx, unsigned int set_lnum);

    double measure() const;    ///< Calculate metrics of the side
    arma::vec3 centre() const; ///< Centre of side
    arma::vec3 normal() const; ///< Vector of (generalized) normal
    double diameter() const;   ///< Calculate the side diameter.

    /// Returns dimension of the side, that is dimension of the element minus one.
    unsigned int dim() const;

    /// Returns true for all sides either on boundary or connected to vb neigboring.
    bool is_external() const;

    /// Returns true for side on the boundary.
    bool is_boundary() const;

    /// Returns node for given local index @p i on the side.
    NodeAccessor<3> node(unsigned int i) const;

    /// Returns iterator to the element of the side.
    ElementAccessor<3> element() const;

    /// Returns global index of the edge connected to the side.
    unsigned int edge_idx() const;

    /// Returns pointer to the edge connected to the side.
    Edge edge() const;

    /** 
     * Returns boundary condition prescribed on the side.
     * Fails on assert if side if not on boundary and no BC is prescribed.
     */ 
    Boundary cond() const;

    /// Returns global index of the prescribed boundary condition.
    unsigned int cond_idx() const;



    /// Returns number of nodes of the side.
    unsigned int n_nodes() const
    { return dim()+1; }

    /// Returns pointer to the mesh.
    const Mesh * mesh() const
    { return this->mesh_; }

    /// Returns local index of the side on the element.
    unsigned int side_idx() const
    { return side_idx_; }

    /// Returns index of element in Mesh::element_vec_.
    unsigned int elem_idx() const
    { return elem_idx_; }

    /// Returns true if the side has assigned element.
    bool is_valid() const
    { return mesh_!= nullptr; }

    /// Iterate over local sides of the element.
    void inc() {
        ASSERT_DBG(is_valid()).error("Do not call inc() for invalid accessor!");
        side_idx_++;
    }

    bool operator==(const Side &other) const {
        return (mesh_ == other.mesh_ ) && ( elem_idx_ == other.elem_idx_ )
        		&& ( side_idx_ == other.side_idx_ );
    }

    bool operator!=(const Side &other) const {
        return !( *this == other);
    }

    /// This is necessary by current DofHandler, should change this
    //void *make_ptr() const;
private:

    arma::vec3 normal_point() const;
    arma::vec3 normal_line() const;
    arma::vec3 normal_triangle() const;

    // Topology of the mesh

    const Mesh * mesh_;     ///< Pointer to Mesh to which belonged
    unsigned int elem_idx_; ///< Index of element in Mesh::element_vec_
    unsigned int side_idx_; ///< Local # of side in element  (to remove it, we heve to remove calc_side_rhs)

};


/*
 * Iterator to a side.
 */
class SideIter {
public:
    SideIter()
    {}

    SideIter(const Side &side)
    : side_(side)
    {}

    bool operator==(const SideIter &other) {
        return (side_.mesh() == other.side_.mesh() ) && ( side_.elem_idx() == other.side_.elem_idx() )
        		&& ( side_.side_idx() == other.side_.side_idx() );
    }


    bool operator!=(const SideIter &other) {
        return !( *this == other);
    }

    ///  * dereference operator
    const Side & operator *() const
            { return side_; }

    /// -> dereference operator
    const Side * operator ->() const
            { return &side_; }

    /// prefix increment iterate only on local element
    SideIter &operator ++ () {
        side_.inc();
        return (*this);
    }

private:
    Side side_;
};



#include "mesh/accessors_impl.hh"

#endif /* ACCESSORS_HH_ */
