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
 * @file    mh_dofhandler.hh
 * @brief   
 */

#ifndef MH_DOFHANDLER_HH_
#define MH_DOFHANDLER_HH_

#include <ext/alloc_traits.h>                // for __alloc_traits<>::value_...
#include <sys/types.h>                       // for uint
#include <memory>                            // for shared_ptr, __shared_ptr
#include <unordered_map>                     // for unordered_map
#include <vector>                            // for vector
#include <armadillo>
#include "la/distribution.hh"                // for Distribution
#include "mesh/long_idx.hh"                  // for LongIdx
#include "mesh/accessors.hh"                 // for ElementAccessor
#include "mesh/elements.h"                   // for Element::side, Element::dim
#include "mesh/mesh.h"                       // for Mesh
#include "mesh/region.hh"                    // for Region
#include "mesh/side_impl.hh"                 // for Side::edge_idx
#include "mesh/sides.h"                      // for SideIter, Side
#include "fem/dh_cell_accessor.hh"           // for DHCellAccessor

class LocalToGlobalMap;
template <int spacedim> class LocalElementAccessorBase;

using namespace std;

/// temporary solution to provide access to results
/// from DarcyFlowMH independent of mesh
class MH_DofHandler {
public:
    MH_DofHandler();
    ~MH_DofHandler();
    void reinit(Mesh *mesh);

    void prepare_parallel();
    void make_row_numberings();
    void prepare_parallel_bddc();

    void set_solution( double time, double * solution, double precision);

    inline double time_changed() const
        { return time_; }

    unsigned int side_dof(const SideIter side) const;

    /// temporary replacement for DofHandler accessor, flux through given side
    double side_flux(const Side &side) const;

    /// temporary replacement for DofHandler accessor, scalar (pressure) on edge of the side
    double side_scalar(const Side &side) const;

    /// temporary replacement for DofHandler accessor, scalar (pressure) on element
    double element_scalar( ElementAccessor<3> &ele ) const;

    inline double precision() const { return solution_precision; };

//protected:
    vector< vector<unsigned int> > elem_side_to_global;

    Mesh *mesh_;
    LongIdx *el_4_loc;          //< array of idexes of local elements (in ordering matching the optimal global)
    LongIdx *row_4_el;          //< element index to matrix row
    LongIdx *side_id_4_loc;     //< array of ids of local sides
    LongIdx *side_row_4_id;     //< side id to matrix row
    LongIdx *edge_4_loc;        //< array of indexes of local edges
    LongIdx *row_4_edge;        //< edge index to matrix row

    // parallel
    Distribution *edge_ds;          //< optimal distribution of edges
    Distribution *el_ds;            //< optimal distribution of elements
    Distribution *side_ds;          //< optimal distribution of elements
    std::shared_ptr<Distribution> rows_ds;          //< final distribution of rows of MH matrix


    /// Maps mesh index of the edge to the edge index in the mesh portion local to the processor.
    /// Temporary solution until we have parallel mesh which should provide such information.
    std::unordered_map<unsigned int, unsigned int> edge_new_local_4_mesh_idx_;

    /// Necessary only for BDDC solver.
    std::shared_ptr<LocalToGlobalMap> global_row_4_sub_row;           //< global dof index for subdomain index


    double * mh_solution;
    double solution_precision;
    double time_;

    friend LocalElementAccessorBase<3>;
};



typedef unsigned int uint;

template <int spacedim>
class LocalElementAccessorBase {
public:

    LocalElementAccessorBase(MH_DofHandler *dh, DHCellAccessor dh_cell)
    : dh(dh), dh_cell_(dh_cell), global_indices_(dh_cell_.dh()->max_elem_dofs()),
	  local_indices_(dh_cell_.dh()->max_elem_dofs())
    {
        n_indices_ = dh_cell_.get_dof_indices(global_indices_);
        dh_cell_.get_loc_dof_indices(local_indices_);
    }

    inline DHCellAccessor dh_cell() const {
        return dh_cell_;
    }

    uint dim() const {
        return dh_cell_.dim();
    }

    uint n_sides() const {
        return element_accessor()->n_sides();
    }

    inline ElementAccessor<3> element_accessor() const {
        return dh_cell_.elm();
    }

    const arma::vec3 centre() const {
        return element_accessor().centre();
    }

    double measure() const {
        return element_accessor().measure();
    }

    Region region() const {
        return element_accessor().region();
    }

    uint ele_global_idx() {
        return element_accessor().idx();
    }

    uint ele_local_idx() const {
        return dh_cell_.local_idx();
    }

    uint ele_row() {
        //return dh->row_4_el[ele_global_idx()];
        return global_indices_[n_indices_/2];
    }

    uint ele_local_row() {
        //return ele_row() - dh->rows_ds->begin(); //  i_loc_el + side_ds->lsize();
        return local_indices_[n_indices_/2];
    }

    uint edge_row(uint i) {
        return global_indices_[(n_indices_+1)/2+i];
    }

    uint edge_local_row( uint i) {
        return local_indices_[(n_indices_+1)/2+i];
    }

    SideIter side(uint i) {
        return element_accessor().side(i);
    }

    uint side_row(uint i) {
        return global_indices_[i];
    }

    uint side_local_row( uint i) {
        return local_indices_[i];
    }

private:
    MH_DofHandler *dh;
    DHCellAccessor dh_cell_;
    std::vector<LongIdx> global_indices_;
    std::vector<LongIdx> local_indices_;
    uint n_indices_;
};

/**
 * This is prototype of further much more complex and general accessor templated by
 * element dimension. In fact we shall need an accessor for every kind of element interaction integral.
 */
template <int spacedim, int dim>
class LocalElementAccessor : public LocalElementAccessorBase<spacedim> {
public:
    LocalElementAccessor(MH_DofHandler & dh, uint loc_ele_idx)
    : LocalElementAccessorBase<spacedim>(dh, loc_ele_idx)
    {}
};



#endif /* MH_DOFHANDLER_HH_ */
