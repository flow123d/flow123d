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
 * @file    integral_acc.hh
 * @brief
 * @author  David Flanderka
 */

#ifndef INTEGRAL_ACC_HH_
#define INTEGRAL_ACC_HH_

#include <memory>
#include <armadillo>
#include "fem/eval_points.hh"
#include "fem/integral_points.hh"
#include "fem/element_cache_map.hh"
#include "mesh/range_wrapper.hh"
#include "mesh/accessors.hh"
#include "fem/dh_cell_accessor.hh"
#include "fem/patch_fe_values.hh"
#include "fem/op_function.hh"
#include "fem/op_accessors_impl.hh"


class BulkIntegral;
class EdgeIntegral;
class CouplingIntegral;
class BoundaryIntegral;
template <unsigned int qdim> class BulkIntegralAcc;
template <unsigned int qdim> class EdgeIntegralAcc;
template <unsigned int qdim> class CouplingIntegralAcc;
template <unsigned int qdim> class BoundaryIntegralAcc;


/**
 * Base integral class holds common data members and methods.
 */
class BaseIntegral {
public:
    /// Default constructor
	BaseIntegral() : eval_points_(nullptr), dim_(0) {}

    /// Constructor of bulk or side subset
	BaseIntegral(std::shared_ptr<EvalPoints> eval_points, unsigned int dim)
	 : eval_points_(eval_points), dim_(dim) {}

    /// Destructor
    virtual ~BaseIntegral();

    /// Getter of eval_points
    std::shared_ptr<EvalPoints> eval_points() const {
        return eval_points_;
    }

    /// Returns dimension.
    unsigned int dim() const {
    	return dim_;
    }
protected:
    /// Pointer to EvalPoints
    std::shared_ptr<EvalPoints> eval_points_;
    /// Dimension of the cell on which points are placed
    unsigned int dim_;
};

/**
 * Temporary class. Parent of all integral accessor.
 *
 * Will be merged with BaseIntegral
 */
template<unsigned int dim>
class FactoryBase
{
protected:
    // Default constructor
    FactoryBase() : patch_fe_values_(nullptr), element_cache_map_(nullptr), quad_(nullptr)
    {}

    // Constructor
    FactoryBase(PatchFEValues<3> *pfev, ElementCacheMap *element_cache_map, std::shared_ptr< FiniteElement<dim> > fe, Quadrature *quad)
    : patch_fe_values_(pfev), element_cache_map_(element_cache_map), fe_(fe), quad_(quad)
    {}

    /// Factory method. Creates operation of given OpType.
    template<class OpType>
    PatchOp<3> *make_patch_op() {
        return patch_fe_values_->get< OpType, dim >(quad_);
    }

    /// Factory method. Same as previous but creates FE operation.
    template<class ValueType, template<unsigned int, class, unsigned int> class OpType, class Domain>
    FeQArray<ValueType> make_qarray(uint component_idx = 0) {
        std::shared_ptr<FiniteElement<dim>> fe_component = patch_fe_values_->fe_comp(fe_, component_idx);
        return FeQArray<ValueType>(patch_fe_values_->template get< OpType<dim, Domain, 3>, dim >(quad_, fe_component));
    }

    PatchFEValues<3> *patch_fe_values_;
    ElementCacheMap *element_cache_map_;
    std::shared_ptr< FiniteElement<dim> > fe_;
    Quadrature *quad_;

    friend class BulkIntegralAcc<dim>;
    friend class EdgeIntegralAcc<dim>;
    friend class CouplingIntegralAcc<dim>;
    friend class BoundaryIntegralAcc<dim>;
};


namespace internal_integrals {

class Base {
public:
    /// Default constructor
	Base() : dim_(0), quad_(nullptr) {}

    /// Constructor of bulk or side subset
	Base(Quadrature *quad, unsigned int dim)
	 : dim_(dim), quad_(quad) {}

    /// Destructor
    virtual ~Base()
    {}

    /// Returns dimension.
    unsigned int dim() const {
    	return dim_;
    }

    /// Returns quadrature.
    Quadrature *quad() const {
    	return quad_;
    }

    /// Comparison operator
    bool operator==(const Base& other) const {
        return (dim_ == other.dim_) && (quad()->size() == other.quad()->size());
    }
protected:
    /// Dimension of the cell on which points are placed
    unsigned int dim_;
    /// Pointer to Quadrature that represents quadrature points of integral.
    Quadrature *quad_;
};

class Bulk : public Base {
public:
    typedef BulkPoint PointType;
    typedef unsigned int MeshItem;

    /// Default constructor
    Bulk() : Base() {}

    /// Constructor of bulk integral- obsolete constructor
    Bulk(Quadrature *quad, unsigned int dim, std::shared_ptr<EvalPoints> eval_points, unsigned int subset_idx)
     : Base(quad, dim) {
        subset_index_ = subset_idx;
        begin_idx_ = eval_points->subset_begin(dim_, subset_index_);
        end_idx_ = eval_points->subset_end(dim_, subset_index_);
    }

    /// Destructor
    ~Bulk()
    {}

    /// Getter of bulk_begin
    uint begin_idx() const {
        return begin_idx_;
    }

private:
    /// Index of data block according to subset in EvalPoints object.
    unsigned int subset_index_;
    uint begin_idx_;
    uint end_idx_;

    friend class ::BulkIntegral;
    friend class ::CouplingIntegral;
    friend class ::BoundaryIntegral;
    template <unsigned int qdim>
    friend class ::BulkIntegralAcc;
    template <unsigned int qdim>
    friend class ::CouplingIntegralAcc;
    template <unsigned int qdim>
    friend class ::BoundaryIntegralAcc;
};

class Edge : public Base {
public:
    typedef EdgePoint PointType;
    typedef DHCellSide MeshItem;

    /// Default constructor
    Edge() : Base() {}

    /// Constructor of edge integral
    Edge(Quadrature *quad, unsigned int dim, std::shared_ptr<EvalPoints> eval_points, unsigned int subset_idx)
     : Base(quad, dim) {
        subset_index_ = subset_idx;

        begin_idx_ = eval_points->subset_begin(dim_, subset_index_);
        uint end_idx = eval_points->subset_end(dim_, subset_index_);
        n_sides_ = dim_ + 1;
        //DebugOut() << "begin: " << begin_idx_ << "end: " << end_idx;
        n_points_per_side_ = (end_idx - begin_idx_) / n_sides_;
        //DebugOut() << "points per side: " << n_points_per_side_;
    }

    /// Destructor
    ~Edge()
    {}

    inline uint side_begin(const DHCellSide &cell_side) const {
        return begin_idx_ + cell_side.side_idx() * n_points_per_side_;
    }

private:
    unsigned int subset_index_;
    uint begin_idx_;

    /// Number of sides (value 0 indicates bulk set)
    unsigned int n_sides_;
    /// Number of points. TODO: pass this to the constructor, avoid extraction from the eval_points
    uint n_points_per_side_;

    friend class ::EdgeIntegral;
    friend class ::CouplingIntegral;
    friend class ::BoundaryIntegral;
    template <unsigned int qdim>
    friend class ::EdgeIntegralAcc;
    template <unsigned int qdim>
    friend class ::CouplingIntegralAcc;
    template <unsigned int qdim>
    friend class ::BoundaryIntegralAcc;
};

}


/**
 * Integral class of bulk points, allows assemblation of volume integrals.
 */
class BulkIntegral : public BaseIntegral {
public:
    /// Default constructor
    BulkIntegral() : BaseIntegral() {}

    /// Constructor of bulk integral
    BulkIntegral(std::shared_ptr<EvalPoints> eval_points, Quadrature *quad, unsigned int dim)
     : BaseIntegral(eval_points, dim)
    {
        switch (dim) {
        case 1:
            internal_bulk_ = eval_points->add_bulk_internal<1>(quad);
            break;
        case 2:
            internal_bulk_ = eval_points->add_bulk_internal<2>(quad);
            break;
        case 3:
            internal_bulk_ = eval_points->add_bulk_internal<3>(quad);
            break;
        default:
            ASSERT_PERMANENT(false);
        }
    }

    /// Destructor
    ~BulkIntegral();

    /// Return index of data block according to subset in EvalPoints object
    inline int get_subset_idx() const {
        return internal_bulk_->subset_index_;
    }


    /// Returns range of bulk local points for appropriate cell accessor - obsolete method
    inline Range< BulkPoint > points(unsigned int element_patch_idx, const ElementCacheMap *elm_cache_map) const {
        auto bgn_it = make_iter<BulkPoint>( BulkPoint(elm_cache_map, element_patch_idx, internal_bulk_->begin_idx_));
        auto end_it = make_iter<BulkPoint>( BulkPoint(elm_cache_map, element_patch_idx, internal_bulk_->end_idx_));
        return Range<BulkPoint>(bgn_it, end_it);
    }

protected:
    /// Internal integral object
    std::shared_ptr<internal_integrals::Bulk> internal_bulk_;
};

/**
 * New Integral accessor class, replace of BulkIntegral, will be merged with BulkIntegral
 *
 * IN DEVELOPMENT
 */
template<unsigned int qdim>
class BulkIntegralAcc : public BulkIntegral {
public:
    /// Default constructor
    BulkIntegralAcc() : BulkIntegral() {}

    /// Constructor of bulk integral
    BulkIntegralAcc(std::shared_ptr<EvalPoints> eval_points, Quadrature *quad, PatchFEValues<3> *pfev, ElementCacheMap *element_cache_map)
     : BulkIntegral(eval_points, quad, qdim), factory_(pfev, element_cache_map, pfev->fe_dim<qdim>(), quad)
    {
        ASSERT_EQ(quad->dim(), qdim);
	}

    /// Destructor
    ~BulkIntegralAcc()
    {}

    /// Returns range of bulk local points for appropriate cell accessor
    inline Range< BulkPoint > points(unsigned int element_patch_idx) const {
        auto bgn_it = make_iter<BulkPoint>( BulkPoint(factory_.element_cache_map_, element_patch_idx, internal_bulk_->begin_idx_));
        auto end_it = make_iter<BulkPoint>( BulkPoint(factory_.element_cache_map_, element_patch_idx, internal_bulk_->end_idx_));
        return Range<BulkPoint>(bgn_it, end_it);
    }

    /**
     * @brief Register the product of Jacobian determinant and the quadrature
     * weight at bulk quadrature points.
     *
     * @param quad Quadrature.
     */
    inline FeQ<Scalar> JxW()
    {
        return FeQ<Scalar>(factory_.template make_patch_op< Op::JxW<qdim, Op::BulkDomain, 3> >());
    }

	/// Create bulk accessor of coords entity
    inline FeQ<Vector> coords()
    {
        return FeQ<Vector>(factory_.template make_patch_op< Op::PtCoords<qdim, Op::BulkDomain, 3> >());
    }

//    inline ElQ<Tensor> jacobian(std::initializer_list<Quadrature *> quad_list)
//    {}

    /// Create bulk accessor of jac determinant entity
    inline ElQ<Scalar> determinant()
    {
        return ElQ<Scalar>( factory_.template make_patch_op< Op::JacDet<qdim, Op::BulkDomain, 3> >() );
    }

    /**
     * @brief Return the value of the @p function_no-th shape function at
     * the @p p bulk quadrature point.
     *
     * @param component_idx Number of the shape function.
     */
    inline FeQArray<Scalar> scalar_shape(uint component_idx = 0)
    {
        return factory_.template make_qarray<Scalar, Op::ScalarShape, Op::BulkDomain>(component_idx);
    }

    /**
     * @brief Return the value of the @p function_no-th vector shape function at
     * the @p p bulk quadrature point.
     *
     * @param component_idx Number of the shape function.
     */
    inline FeQArray<Vector> vector_shape(uint component_idx = 0)
    {
        return factory_.template make_qarray<Vector, Op::DispatchVectorShape, Op::BulkDomain>(component_idx);
    }

//    inline FeQArray<Tensor> tensor_shape(uint component_idx = 0)
//    {}

    /**
     * @brief Return the value of the @p function_no-th gradient shape function at
     * the @p p bulk quadrature point.
     *
     * @param component_idx Number of the shape function.
     */
    inline FeQArray<Vector> grad_scalar_shape(uint component_idx=0)
    {
        return factory_.template make_qarray<Vector, Op::GradScalarShape, Op::BulkDomain>(component_idx);
    }

    /**
     * @brief Return the value of the @p function_no-th gradient vector shape function
     * at the @p p bulk quadrature point.
     *
     * @param component_idx Number of the shape function.
     */
    inline FeQArray<Tensor> grad_vector_shape(uint component_idx=0)
    {
        return factory_.template make_qarray<Tensor, Op::DispatchGradVectorShape, Op::BulkDomain>(component_idx);
    }

    /**
     * @brief Return the value of the @p function_no-th vector symmetric gradient
     * at the @p p bulk quadrature point.
     *
     * @param component_idx Number of the shape function.
     */
    inline FeQArray<Tensor> vector_sym_grad(uint component_idx=0)
    {
        return factory_.template make_qarray<Tensor, Op::VectorSymGrad, Op::BulkDomain>(component_idx);
    }

    /**
     * @brief Return the value of the @p function_no-th vector divergence at
     * the @p p bulk quadrature point.
     *
     * @param component_idx Number of the shape function.
     */
    inline FeQArray<Scalar> vector_divergence(uint component_idx=0)
    {
        return factory_.template make_qarray<Scalar, Op::VectorDivergence, Op::BulkDomain>(component_idx);
    }

private:
    /// Defines interface of operation accessors declaration
    FactoryBase<qdim> factory_;
};

/**
 * Integral class of side points, allows assemblation of element - element fluxes.
 */
class EdgeIntegral : public BaseIntegral {
public:
    /// Default constructor
	EdgeIntegral() : BaseIntegral()
    {
	    ASSERT_PERMANENT(false);
    }

    /// Constructor of edge integral
	EdgeIntegral(std::shared_ptr<EvalPoints> eval_points, Quadrature *quad, unsigned int dim)
	: BaseIntegral(eval_points, dim)
	{
	    switch (dim) {
	    case 1:
	        internal_edge_ = eval_points->add_edge_internal<1>(quad);
	        break;
	    case 2:
	        internal_edge_ = eval_points->add_edge_internal<2>(quad);
	        break;
	    case 3:
	        internal_edge_ = eval_points->add_edge_internal<3>(quad);
	        break;
	    default:
	        ASSERT(false).error("Should not happen!\n");
	    }
	}

    /// Destructor
    ~EdgeIntegral();

    /// Getter of n_sides
    inline unsigned int n_sides() const {
        return internal_edge_->n_sides_;
    }

    /// Return index of data block according to subset in EvalPoints object
    inline int get_subset_idx() const {
        return internal_edge_->subset_index_;
    }

    inline uint side_begin(const DHCellSide &cell_side) const {
        return internal_edge_->side_begin(cell_side);
    }

    /// Returns range of side local points for appropriate cell side accessor - obsolete method
    inline Range< EdgePoint > points(const DHCellSide &cell_side, const ElementCacheMap *elm_cache_map) const {
        ASSERT_EQ(cell_side.dim(), dim_);

        //DebugOut() << "points per side: " << internal_edge_->n_points_per_side_;
        uint element_patch_idx = elm_cache_map->position_in_cache(cell_side.element().idx());
        uint begin_idx = internal_edge_->side_begin(cell_side);
        auto bgn_it = make_iter<EdgePoint>( EdgePoint(
                BulkPoint(elm_cache_map, element_patch_idx, 0), internal_edge_, begin_idx));
        auto end_it = make_iter<EdgePoint>( EdgePoint(
                BulkPoint(elm_cache_map, element_patch_idx, internal_edge_->n_points_per_side_), internal_edge_, begin_idx));
        return Range<EdgePoint>(bgn_it, end_it);
    }


protected:
    /// Internal integral object
    std::shared_ptr<internal_integrals::Edge> internal_edge_;

    friend class EvalPoints;
    friend class EdgePoint;
    friend class CouplingPoint;
    friend class BoundaryPoint;
    friend class CouplingIntegral;
    friend class BoundaryIntegral;
    template<unsigned int qdim>
    friend class CouplingIntegralAcc;
    template<unsigned int qdim>
    friend class BoundaryIntegralAcc;
};

/**
 * New Integral accessor class, replace of EdgeIntegral, will be merged with EdgeIntegral
 *
 * IN DEVELOPMENT
 */
template<unsigned int qdim>
class EdgeIntegralAcc : public EdgeIntegral {
public:
    /// Default constructor
    EdgeIntegralAcc() : EdgeIntegral() {}

    /// Constructor of edge integral
    EdgeIntegralAcc(std::shared_ptr<EvalPoints> eval_points, Quadrature *quad, PatchFEValues<3> *pfev, ElementCacheMap *element_cache_map)
    : EdgeIntegral(eval_points, quad, qdim), factory_(pfev, element_cache_map, pfev->fe_dim<qdim>(), quad)
    {
        ASSERT_EQ(quad->dim()+1, qdim);
    }


    /// Destructor
    ~EdgeIntegralAcc()
    {}

    /// Returns range of side local points for appropriate cell side accessor
    inline Range< EdgePoint > points(const DHCellSide &cell_side) const {
        ASSERT_EQ(cell_side.dim(), dim_);

        //DebugOut() << "points per side: " << internal_edge_->n_points_per_side_;
        uint element_patch_idx = factory_.element_cache_map_->position_in_cache(cell_side.element().idx());
        uint begin_idx = internal_edge_->side_begin(cell_side);
        auto bgn_it = make_iter<EdgePoint>( EdgePoint(
                BulkPoint(factory_.element_cache_map_, element_patch_idx, 0), this->internal_edge_, begin_idx));
        auto end_it = make_iter<EdgePoint>( EdgePoint(
                BulkPoint(factory_.element_cache_map_, element_patch_idx, internal_edge_->n_points_per_side_), this->internal_edge_, begin_idx));
        return Range<EdgePoint>(bgn_it, end_it);
    }

    /// Same as BulkValues::JxW but register at side quadrature points.
    inline FeQ<Scalar> JxW()
    {
        return FeQ<Scalar>(factory_.template make_patch_op< Op::JxW<qdim, Op::SideDomain, 3> >());
    }

    /**
     * @brief Register the normal vector to a side at side quadrature points.
     *
     * @param quad Quadrature.
     */
	inline ElQ<Vector> normal_vector()
	{
        return ElQ<Vector>(factory_.template make_patch_op< Op::NormalVec<qdim, 3> >());
	}

	/// Create side accessor of coords entity
    inline FeQ<Vector> coords()
    {
        return FeQ<Vector>(factory_.template make_patch_op< Op::PtCoords<qdim, Op::SideDomain, 3> >());
    }

    /// Create bulk accessor of jac determinant entity
    inline ElQ<Scalar> determinant()
    {
        return ElQ<Scalar>(factory_.template make_patch_op< Op::JacDet<qdim, Op::SideDomain, 3> >());
    }

    /// Same as BulkValues::scalar_shape but register at side quadrature points.
    inline FeQArray<Scalar> scalar_shape(uint component_idx = 0)
    {
        return factory_.template make_qarray<Scalar, Op::ScalarShape, Op::SideDomain>(component_idx);
    }

    /// Same as BulkValues::vector_shape but register at side quadrature points.
    inline FeQArray<Vector> vector_shape(uint component_idx = 0)
    {
        return factory_.template make_qarray<Vector, Op::DispatchVectorShape, Op::SideDomain>(component_idx);
    }

    /// Same as BulkValues::grad_scalar_shape but register at side quadrature points.
    inline FeQArray<Vector> grad_scalar_shape(uint component_idx=0)
    {
        return factory_.template make_qarray<Vector, Op::GradScalarShape, Op::SideDomain>(component_idx);
    }

    /**
     * @brief Return the value of the @p function_no-th gradient vector shape function
     * at the @p p bulk quadrature point.
     *
     * @param component_idx Number of the shape function.
     */
    inline FeQArray<Tensor> grad_vector_shape(uint component_idx=0)
    {
        return factory_.template make_qarray<Tensor, Op::DispatchGradVectorShape, Op::SideDomain>(component_idx);
    }

    /**
     * @brief Return the value of the @p function_no-th vector symmetric gradient
     * at the @p p side quadrature point.
     *
     * @param component_idx Number of the shape function.
     */
    inline FeQArray<Tensor> vector_sym_grad(uint component_idx=0)
    {
        return factory_.template make_qarray<Tensor, Op::VectorSymGrad, Op::SideDomain>(component_idx);
    }

    /**
     * @brief Return the value of the @p function_no-th vector divergence at
     * the @p p side quadrature point.
     *
     * @param component_idx Number of the shape function.
     */
    inline FeQArray<Scalar> vector_divergence(uint component_idx=0)
    {
        return factory_.template make_qarray<Scalar, Op::VectorDivergence, Op::SideDomain>(component_idx);
    }

private:
    /// Defines interface of operation accessors declaration
    FactoryBase<qdim> factory_;

    friend class EvalPoints;
    friend class EdgePoint;
    friend class CouplingPoint;
    friend class BoundaryPoint;
    template <unsigned int quaddim>
    friend class CouplingIntegralAcc;
    template <unsigned int quaddim>
    friend class BoundaryIntegralAcc;
};


/**
 * Integral class of neighbour points, allows assemblation of element - side fluxes.
 *
 * Dimension corresponds with element of higher dim.
 */
class CouplingIntegral : public BaseIntegral {
public:
    /// Default constructor
	CouplingIntegral() : BaseIntegral() {}

    /// Constructor of ngh integral
	CouplingIntegral(std::shared_ptr<EvalPoints> eval_points, Quadrature *quad, unsigned int dim)
	 : BaseIntegral(eval_points, dim)
	{
	    switch (dim) {
	    case 2:
	        internal_bulk_ = eval_points->add_bulk_internal<1>(quad);
	        internal_edge_ = eval_points->add_edge_internal<2>(quad);
	        break;
	    case 3:
	        internal_bulk_ = eval_points->add_bulk_internal<2>(quad);
	        internal_edge_ = eval_points->add_edge_internal<3>(quad);
	        break;
	    default:
	        ASSERT(false)(dim).error("Should not happen!\n");
	    }
	}

    /// Destructor
    ~CouplingIntegral();

    /// Return index of data block according to subset of higher dim in EvalPoints object
    inline int get_subset_high_idx() const {
        return internal_edge_->subset_index_;
    }

    /// Return index of data block according to subset of lower dim in EvalPoints object
    inline int get_subset_low_idx() const {
        return internal_bulk_->subset_index_;
    }

    inline uint bulk_begin() const {
        return eval_points_->subset_begin(dim_-1, internal_bulk_->subset_index_);
    }

    /// Returns range of side local points for appropriate cell side accessor - obsolete method
    inline Range< CouplingPoint > points(const DHCellSide &cell_side, const ElementCacheMap *elm_cache_map) const {
        ASSERT_EQ(cell_side.dim(), dim_);
        uint element_patch_idx = elm_cache_map->position_in_cache(cell_side.element().idx());
        uint side_begin = internal_edge_->side_begin(cell_side);
        auto bgn_it = make_iter<CouplingPoint>( CouplingPoint(
                BulkPoint(elm_cache_map, element_patch_idx, 0), internal_bulk_, side_begin) );
        auto end_it = make_iter<CouplingPoint>( CouplingPoint(
                BulkPoint(elm_cache_map, element_patch_idx, internal_edge_->n_points_per_side_), internal_bulk_, side_begin) );;
        return Range<CouplingPoint>(bgn_it, end_it);
    }

protected:
    /// Integral according to bulk subset part (element of lower dim) in EvalPoints object.
    std::shared_ptr<internal_integrals::Bulk> internal_bulk_;
    /// Integral according to side subset part (element of higher dim) in EvalPoints object.
    std::shared_ptr<internal_integrals::Edge> internal_edge_;

    friend class CouplingPoint;
};

template<unsigned int qdim>
class CouplingIntegralAcc : public CouplingIntegral {
public:
    /// Default constructor
    CouplingIntegralAcc() : CouplingIntegral() {}

    /// Constructor of ngh integral
    CouplingIntegralAcc(std::shared_ptr<EvalPoints> eval_points, Quadrature *quad, PatchFEValues<3> *pfev, ElementCacheMap *element_cache_map)
     : CouplingIntegral(eval_points, quad, qdim),
	   factory_(pfev, element_cache_map, pfev->fe_dim<qdim>(), quad)
    {
        fe_low_ = pfev->fe_dim<qdim-1>();
    }

    /// Destructor
    ~CouplingIntegralAcc()
    {
    	internal_bulk_.reset();
    	internal_edge_.reset();
    }

    /// Return index of data block according to subset of higher dim in EvalPoints object
    inline int get_subset_high_idx() const {
        return internal_edge_->subset_index_;
    }

    /// Return index of data block according to subset of lower dim in EvalPoints object
    inline int get_subset_low_idx() const {
        return internal_bulk_->subset_index_;
    }

    inline uint bulk_begin() const {
        return internal_bulk_->begin_idx_;
    }

    /// Returns range of side local points for appropriate cell side accessor
    inline Range< CouplingPoint > points(const DHCellSide &cell_side) const {
        ASSERT_EQ(cell_side.dim(), dim_);

        uint element_patch_idx = factory_.element_cache_map_->position_in_cache(cell_side.element().idx());
        uint side_begin = internal_edge_->side_begin(cell_side);
        auto bgn_it = make_iter<CouplingPoint>( CouplingPoint(
                BulkPoint(factory_.element_cache_map_, element_patch_idx, 0), internal_bulk_, side_begin) );
        auto end_it = make_iter<CouplingPoint>( CouplingPoint(
                BulkPoint(factory_.element_cache_map_, element_patch_idx, internal_edge_->n_points_per_side_), internal_bulk_, side_begin) );;
        return Range<CouplingPoint>(bgn_it, end_it);
    }

    /// Factory method. Same as previous but creates FE operation.
    template<class ValueType, template<unsigned int, class, unsigned int> class OpType>
    FeQJoin<ValueType> make_qjoin(uint component_idx = 0) {
        // element of lower dim (bulk points)
        auto fe_component_low = factory_.patch_fe_values_->fe_comp(fe_low_, component_idx);
        auto *low_dim_op = factory_.patch_fe_values_->template get< OpType<qdim-1, Op::BulkDomain, 3>, qdim-1 >(factory_.quad_, fe_component_low);
        auto *low_dim_zero_op = factory_.patch_fe_values_->template get< Op::OpZero<qdim-1, Op::BulkDomain, 3>, qdim-1 >(factory_.quad_, fe_component_low);

    	// element of higher dim (side points)
        auto fe_component_high = factory_.patch_fe_values_->fe_comp(factory_.fe_, component_idx);
        auto *high_dim_op = factory_.patch_fe_values_->template get< OpType<qdim, Op::SideDomain, 3>, qdim >(factory_.quad_, fe_component_high);
        auto *high_dim_zero_op = factory_.patch_fe_values_->template get< Op::OpZero<qdim, Op::SideDomain, 3>, qdim >(factory_.quad_, fe_component_high);

        ASSERT_EQ(fe_component_high->fe_type(), fe_component_low->fe_type()).error("Type of FiniteElement of low and high element must be same!\n");
        return FeQJoin<ValueType>(low_dim_op, high_dim_op, low_dim_zero_op, high_dim_zero_op);
    }

    /// Same as BulkValues::JxW but register at side quadrature points.
    inline FeQ<Scalar> JxW()
    {
        return FeQ<Scalar>(factory_.template make_patch_op< Op::JxW<qdim, Op::SideDomain, 3> >());
    }

    /**
     * @brief Register the normal vector to a side at side quadrature points.
     *
     * @param quad Quadrature.
     */
	inline ElQ<Vector> normal_vector()
	{
        return ElQ<Vector>(factory_.template make_patch_op< Op::NormalVec<qdim, 3> >());
	}

    inline FeQArray<Vector> vector_shape(uint component_idx = 0)
    {
        return factory_.template make_qarray<Vector, Op::DispatchVectorShape, Op::SideDomain>(component_idx);
    }

    inline FeQJoin<Scalar> scalar_join_shape(uint component_idx = 0)
    {
        return this->template make_qjoin<Scalar, Op::ScalarShape>(component_idx);
    }

    inline FeQJoin<Vector> vector_join_shape(uint component_idx = 0)
    {
        return this->template make_qjoin<Vector, Op::DispatchVectorShape>(component_idx);
    }

    inline FeQJoin<Tensor> gradient_vector_join_shape(uint component_idx = 0)
    {
        return this->template make_qjoin<Tensor, Op::DispatchGradVectorShape>(component_idx);
    }


private:
    /// Defines interface of operation accessors declaration
    FactoryBase<qdim> factory_;

    /// Holds FiniteEementt object of lower dimension
    std::shared_ptr< FiniteElement<qdim-1> > fe_low_;

    friend class CouplingPoint;
};

/// Template specialization of previous class
template<>
class CouplingIntegralAcc<1> : public CouplingIntegral {
public:
    /// Default constructor
    CouplingIntegralAcc() : CouplingIntegral() {}

    /// Constructor of ngh integral
    CouplingIntegralAcc(std::shared_ptr<EvalPoints> eval_points, Quadrature *quad, PatchFEValues<3> *pfev, ElementCacheMap *element_cache_map)
     : CouplingIntegral(),
	   factory_(pfev, element_cache_map, pfev->fe_dim<1>(), quad)
    {
        ASSERT_EQ(quad->dim(), 0);
        this->eval_points_ = eval_points;
        this->dim_ = 1;
    }

    /// Destructor
    ~CouplingIntegralAcc()
    {}

    /// Returns empty point range
    inline Range< CouplingPoint > points(const DHCellSide &cell_side) const {
        ASSERT_EQ(cell_side.dim(), dim_);
        auto iter = make_iter<CouplingPoint>( CouplingPoint() );
        return Range<CouplingPoint>(iter, iter);
    }

    /// Define empty operations
    inline FeQ<Scalar> JxW()
    {
        return FeQ<Scalar>();
    }

    inline ElQ<Vector> normal_vector()
    {
        return ElQ<Vector>();
    }

    inline FeQArray<Vector> vector_shape(FMT_UNUSED uint component_idx = 0)
    {
        return FeQArray<Vector>();
    }

    inline FeQJoin<Scalar> scalar_join_shape(FMT_UNUSED uint component_idx = 0)
    {
        return FeQJoin<Scalar>();
    }

    inline FeQJoin<Vector> vector_join_shape(FMT_UNUSED uint component_idx = 0)
    {
        return FeQJoin<Vector>();
    }

    inline FeQJoin<Tensor> gradient_vector_join_shape(FMT_UNUSED uint component_idx = 0)
    {
        return FeQJoin<Tensor>();
    }

private:
    /// Defines interface of operation accessors declaration
    FactoryBase<1> factory_;
};


/**
 * Integral class of boundary points, allows assemblation of fluxes between sides and neighbouring boundary elements.
 */
class BoundaryIntegral : public BaseIntegral {
public:
    /// Default constructor
    BoundaryIntegral() : BaseIntegral() {}

    /// Constructor of bulk subset
    BoundaryIntegral(std::shared_ptr<EvalPoints> eval_points, Quadrature *quad, unsigned int dim);

    /// Destructor
    ~BoundaryIntegral();

    /// Return index of data block according to subset of higher dim in EvalPoints object
    inline int get_subset_high_idx() const {
        return internal_edge_->subset_index_;
    }

    /// Return index of data block according to subset of lower dim (boundary) in EvalPoints object
    inline int get_subset_low_idx() const {
        return internal_bulk_->subset_index_;
    }

    inline uint bulk_begin() const {
      //  DebugOut().fmt("edge_begin: {} bdr_begin: {}",
      //          internal_bulk_->get_begin_idx(),
      //          internal_bulk_->get_begin_idx());
        return internal_bulk_->begin_idx_;
    }

    /// Returns range of bulk local points for appropriate cell accessor - obsolete method
    inline Range< BoundaryPoint > points(const DHCellSide &cell_side, const ElementCacheMap *elm_cache_map) const {
        ASSERT_EQ(cell_side.dim(), dim_);
        uint element_patch_idx = elm_cache_map->position_in_cache(cell_side.element().idx());
        uint side_begin = internal_edge_->side_begin(cell_side);
        auto bgn_it = make_iter<BoundaryPoint>( BoundaryPoint(
                BulkPoint(elm_cache_map, element_patch_idx, 0), internal_bulk_, side_begin) );
        auto end_it = make_iter<BoundaryPoint>( BoundaryPoint(
                BulkPoint(elm_cache_map, element_patch_idx, internal_edge_->n_points_per_side_), internal_bulk_, side_begin) );;
        return Range<BoundaryPoint>(bgn_it, end_it);
    }

protected:
    /// Integral according to kower dim (boundary) element subset part in EvalPoints object.
    std::shared_ptr<internal_integrals::Bulk> internal_bulk_;
    /// Integral according to higher dim (bulk) element subset part in EvalPoints object.
    std::shared_ptr<internal_integrals::Edge> internal_edge_;

    friend class BoundaryPoint;
};


template<unsigned int qdim>
class BoundaryIntegralAcc : public BoundaryIntegral {
public:
    /// Default constructor
    BoundaryIntegralAcc() : BoundaryIntegral() {}

    BoundaryIntegralAcc(std::shared_ptr<EvalPoints> eval_points, Quadrature *quad, PatchFEValues<3> *pfev, ElementCacheMap *element_cache_map)
    : BoundaryIntegral(eval_points, quad, qdim), factory_(pfev, element_cache_map, pfev->fe_dim<qdim>(), quad)
    {
        ASSERT_EQ(quad->dim()+1, qdim);
    }

    /// Destructor
    ~BoundaryIntegralAcc()
    {}

    /// Returns range of bulk local points for appropriate cell accessor
    inline Range< BoundaryPoint > points(const DHCellSide &cell_side) const {
        ASSERT_EQ(cell_side.dim(), dim_);

        uint element_patch_idx = factory_.element_cache_map_->position_in_cache(cell_side.element().idx());
        uint side_begin = internal_edge_->side_begin(cell_side);
        auto bgn_it = make_iter<BoundaryPoint>( BoundaryPoint(
                BulkPoint(factory_.element_cache_map_, element_patch_idx, 0), internal_bulk_, side_begin) );
        auto end_it = make_iter<BoundaryPoint>( BoundaryPoint(
                BulkPoint(factory_.element_cache_map_, element_patch_idx, internal_edge_->n_points_per_side_), internal_bulk_, side_begin) );;
        return Range<BoundaryPoint>(bgn_it, end_it);
    }

    /// Same as BulkValues::JxW but register at side quadrature points.
    inline FeQ<Scalar> JxW()
    {
        return FeQ<Scalar>(factory_.template make_patch_op< Op::JxW<qdim, Op::SideDomain, 3> >());
    }

    /**
     * @brief Register the normal vector to a side at side quadrature points.
     *
     * @param quad Quadrature.
     */
	inline ElQ<Vector> normal_vector()
	{
        return ElQ<Vector>(factory_.template make_patch_op< Op::NormalVec<qdim, 3> >());
	}

	/// Create side accessor of coords entity
    inline FeQ<Vector> coords()
    {
        return FeQ<Vector>(factory_.template make_patch_op< Op::PtCoords<qdim, Op::SideDomain, 3> >());
    }

    /// Create bulk accessor of jac determinant entity
    inline ElQ<Scalar> determinant()
    {
        return ElQ<Scalar>(factory_.template make_patch_op< Op::JacDet<qdim, Op::SideDomain, 3> >());
    }

    /// Same as BulkValues::scalar_shape but register at side quadrature points.
    inline FeQArray<Scalar> scalar_shape(uint component_idx = 0)
    {
        return factory_.template make_qarray<Scalar, Op::ScalarShape, Op::SideDomain>(component_idx);
    }

    /// Same as BulkValues::vector_shape but register at side quadrature points.
    inline FeQArray<Vector> vector_shape(uint component_idx = 0)
    {
        return factory_.template make_qarray<Vector, Op::DispatchVectorShape, Op::SideDomain>(component_idx);
    }

    /// Same as BulkValues::grad_scalar_shape but register at side quadrature points.
    inline FeQArray<Vector> grad_scalar_shape(uint component_idx=0)
    {
        return factory_.template make_qarray<Vector, Op::GradScalarShape, Op::SideDomain>(component_idx);
    }

    /**
     * @brief Return the value of the @p function_no-th gradient vector shape function
     * at the @p p bulk quadrature point.
     *
     * @param component_idx Number of the shape function.
     */
    inline FeQArray<Tensor> grad_vector_shape(uint component_idx=0)
    {
        return factory_.template make_qarray<Tensor, Op::DispatchGradVectorShape, Op::SideDomain>(component_idx);
    }

    /**
     * @brief Return the value of the @p function_no-th vector symmetric gradient
     * at the @p p side quadrature point.
     *
     * @param component_idx Number of the shape function.
     */
    inline FeQArray<Tensor> vector_sym_grad(uint component_idx=0)
    {
        return factory_.template make_qarray<Tensor, Op::VectorSymGrad, Op::SideDomain>(component_idx);
    }

    /**
     * @brief Return the value of the @p function_no-th vector divergence at
     * the @p p side quadrature point.
     *
     * @param component_idx Number of the shape function.
     */
    inline FeQArray<Scalar> vector_divergence(uint component_idx=0)
    {
        return factory_.template make_qarray<Scalar, Op::VectorDivergence, Op::SideDomain>(component_idx);
    }

private:
    /// Defines interface of operation accessors declaration
    FactoryBase<qdim> factory_;

    friend class BoundaryPoint;
};


#endif /* INTEGRAL_ACC_HH_ */
