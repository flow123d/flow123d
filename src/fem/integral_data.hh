/*!
 *
? * Copyright (C) 2015 Technical University of Liberec.  All rights reserved.
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
 * @file    integral_data.hh
 * @brief
 * @author  David Flanderka
 */

#ifndef INTEGRAL_DATA_HH_
#define INTEGRAL_DATA_HH_

#include <vector>
#include <memory>
#include <armadillo>
#include <unordered_set>
#include <unordered_map>
#include <boost/functional/hash.hpp>      // for boost::hash_value
#include "fem/dh_cell_accessor.hh"
#include "tools/revertable_list.hh"
#include "mesh/range_wrapper.hh"

template <unsigned int dim> class BulkIntegralAcc;
template <unsigned int dim> class EdgeIntegralAcc;
template <unsigned int dim> class CouplingIntegralAcc;
template <unsigned int dim> class BoundaryIntegralAcc;
template <int spacedim> class ElementAccessor;



/**
 * Helper structure holds data of cell (bulk) integral
 *
 * Data is specified by cell and subset index in EvalPoint object
 */
struct BulkIntegralData {
	/// Default constructor
    BulkIntegralData() {}

    /// Constructor with data mebers initialization
    BulkIntegralData(DHCellAccessor dhcell)
    : cell(dhcell) {}

    /// Copy constructor
    BulkIntegralData(const BulkIntegralData &other)
    : cell(other.cell) {}

    DHCellAccessor cell;          ///< Specified cell (element)
};

/**
 * Helper structure holds data of edge integral
 *
 * Data is specified by side and subset index in EvalPoint object
 */
struct EdgeIntegralData {
	/// Default constructor
	EdgeIntegralData()
	: edge_side_range(make_iter<DHEdgeSide, DHCellSide>( DHEdgeSide() ), make_iter<DHEdgeSide, DHCellSide>( DHEdgeSide() )) {}

    /// Copy constructor
	EdgeIntegralData(const EdgeIntegralData &other)
    : edge_side_range(other.edge_side_range) {}

    /// Constructor with data mebers initialization
	EdgeIntegralData(RangeConvert<DHEdgeSide, DHCellSide> range)
    : edge_side_range(range) {}

	RangeConvert<DHEdgeSide, DHCellSide> edge_side_range;   ///< Specified cell side (element)
};

/**
 * Helper structure holds data of neighbour (coupling) integral
 *
 * Data is specified by cell, side and their subset indices in EvalPoint object
 */
struct CouplingIntegralData {
	/// Default constructor
   	CouplingIntegralData() {}

    /// Constructor with data mebers initialization
   	CouplingIntegralData(DHCellAccessor dhcell, DHCellSide dhside)
    : cell(dhcell), side(dhside) {}

    /// Copy constructor
   	CouplingIntegralData(const CouplingIntegralData &other)
    : cell(other.cell), side(other.side) {}

    DHCellAccessor cell;
    DHCellSide side;                   ///< Specified cell side (higher dim element)
};

/**
 * Helper structure holds data of boundary integral
 *
 * Data is specified by side and subset indices of side and appropriate boundary element in EvalPoint object
 */
struct BoundaryIntegralData {
	/// Default constructor
	BoundaryIntegralData() {}

    /// Constructor with data mebers initialization
	BoundaryIntegralData(DHCellSide dhside)
    : side(dhside) {}

    /// Copy constructor
	BoundaryIntegralData(const BoundaryIntegralData &other)
    : side(other.side) {}

	// We don't need hold ElementAccessor of boundary element, side.cond().element_accessor() provides it.
    DHCellSide side;                   ///< Specified cell side (bulk element)
};




/// Define Integral Tuple hash function - helper struct of IntegralPtrMap
struct IntegralTplHash {
    std::size_t operator()(std::tuple<uint, uint> tpl) const {
        return boost::hash_value( tpl );
    }

    /// Create tuple from dimennsion and size of Quadrature
    static std::tuple<uint, uint> integral_tuple(uint dim, uint quad_size) {
        return std::make_tuple(dim, quad_size);
    }

};

/// Alias for unordered_map of shared_ptr<Integral> with custom hash
template<typename Integral>
using IntegralPtrMap = std::unordered_map<std::tuple<uint, uint>, std::shared_ptr<Integral>, IntegralTplHash>;


/// Set of integral of given dimension necessary in assemblation
template<unsigned int dim>
struct DimIntegrals {
public:
	IntegralPtrMap<BulkIntegralAcc<dim>> bulk_;            ///< Bulk integrals of elements
	IntegralPtrMap<EdgeIntegralAcc<dim>> edge_;            ///< Edge integrals between elements of same dimensions
	IntegralPtrMap<CouplingIntegralAcc<dim>> coupling_;    ///< Coupling integrals between elements of dimensions dim and dim-1
	IntegralPtrMap<BoundaryIntegralAcc<dim>> boundary_;    ///< Boundary integrals betwwen side and boundary element of dim-1


    /// Finalize temporary part of integral data.
    inline void make_permanent()
    {
        for (auto &it : bulk_) it.second->patch_data().make_permanent();
        for (auto &it : edge_) it.second->patch_data().make_permanent();
        for (auto &it : coupling_) it.second->patch_data().make_permanent();
        for (auto &it : boundary_) it.second->patch_data().make_permanent();
    }

    /// Erase temporary part of integral data.
    inline void revert_temporary()
    {
        for (auto &it : bulk_) it.second->patch_data().revert_temporary();
        for (auto &it : edge_) it.second->patch_data().revert_temporary();
        for (auto &it : coupling_) it.second->patch_data().revert_temporary();
        for (auto &it : boundary_) it.second->patch_data().revert_temporary();
    }

    /// Clear list of integral data.
    inline void reset()
    {
        for (auto &it : bulk_) it.second->patch_data().reset();
        for (auto &it : edge_) it.second->patch_data().reset();
        for (auto &it : coupling_) it.second->patch_data().reset();
        for (auto &it : boundary_) it.second->patch_data().reset();
    }

    /// Return number of cells on patch
    inline unsigned int n_patch_cells() const {
        if (bulk_.size() > 0) return bulk_.begin()->second->patch_data().permanent_size();
        else return 0;
    }

    /// Return number of edges on patch
    inline unsigned int n_patch_edges() const {
        if (edge_.size() > 0) return edge_.begin()->second->patch_data().permanent_size();
        else return 0;
    }

    /// Return number of neighbours on patch
    inline unsigned int n_patch_neighbours() const {
        if (coupling_.size() > 0) return coupling_.begin()->second->patch_data().permanent_size();
        else return 0;
    }

    /// Return number of boundaries on patch
    inline unsigned int n_patch_boundaries() const {
        if (boundary_.size() > 0) return boundary_.begin()->second->patch_data().permanent_size();
        else return 0;
    }

};




/// Define Integral Tuple hash function - helper struct of OperationMap
struct OperationTplHash {
    std::size_t operator()(std::tuple<std::string, uint> tpl) const {
        return boost::hash_value( tpl );
    }

    /// Create tuple from typeid(Operation).name and size of Quadrature
    static std::tuple<std::string, uint> op_tuple(std::string op_type, uint quad_size) {
        return std::make_tuple(op_type, quad_size);
    }

};

/// Alias for unordered_map of Operation pointer with custom hash
template<typename Operation>
using OperationMap = std::unordered_map<std::tuple<std::string, uint>, Operation *, OperationTplHash>;



#endif /* INTEGRAL_DATA_HH_ */
