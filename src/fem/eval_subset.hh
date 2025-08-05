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
 * @file    eval_subset.hh
 * @brief
 * @author  David Flanderka
 */

#ifndef EVAL_SUBSET_HH_
#define EVAL_SUBSET_HH_

#include <memory>
#include <armadillo>
#include "fem/eval_points.hh"
#include "fem/eval_subset_points.hh"
#include "fem/element_cache_map.hh"
#include "mesh/range_wrapper.hh"
#include "mesh/accessors.hh"
#include "fem/dh_cell_accessor.hh"
#include "fem/patch_fe_values.hh"
#include "fem/op_function.hh"
#include "fem/op_accessors_impl.hh"


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
public:
    /// Getter of quadrature
    const Quadrature *get_quad() const {
        return quad_;
    }

protected:
    // Default constructor
    FactoryBase() : patch_fe_values_(nullptr), quad_(nullptr)
    {}

    // Constructor
    FactoryBase(PatchFEValues<3> *pfev, Quadrature *quad) : patch_fe_values_(pfev), quad_(quad)
    {}

    /// Factory method. Creates operation of given OpType.
    template<class OpType>
    PatchOp<3> *make_patch_op() {
        return patch_fe_values_->get< OpType, dim >();
    }

    /// Factory method. Same as previous but creates FE operation.
    template<class ValueType, template<unsigned int, class, unsigned int> class OpType, class Domain>
    FeQArray<ValueType> make_qarray(uint component_idx = 0) {
        std::shared_ptr<FiniteElement<dim>> fe_component = patch_fe_values_->fe_comp(fe_, component_idx);
        return FeQArray<ValueType>(patch_fe_values_->template get< OpType<dim, Domain, 3>, dim >(fe_component));
    }

    PatchFEValues<3> *patch_fe_values_;
    std::shared_ptr< FiniteElement<dim> > fe_;
    Quadrature *quad_;
};


/**
 * Integral class of bulk points, allows assemblation of volume integrals.
 */
class BulkIntegral : public BaseIntegral, public std::enable_shared_from_this<BulkIntegral> {
public:
    /// Default constructor
	BulkIntegral() : BaseIntegral() {}

    /// Constructor of bulk integral
	BulkIntegral(std::shared_ptr<EvalPoints> eval_points, unsigned int dim, uint i_subset)
	 : BaseIntegral(eval_points, dim), subset_index_(i_subset)
	{
	    begin_idx_ = eval_points_->subset_begin(dim_, subset_index_);
	    end_idx_ = eval_points_->subset_end(dim_, subset_index_);
	}

    /// Destructor
    ~BulkIntegral();

    /// Return index of data block according to subset in EvalPoints object
    inline int get_subset_idx() const {
        return subset_index_;
    }


    /// Returns range of bulk local points for appropriate cell accessor
    inline Range< BulkPoint > points(unsigned int element_patch_idx, const ElementCacheMap *elm_cache_map) const {
        auto bgn_it = make_iter<BulkPoint>( BulkPoint(elm_cache_map, element_patch_idx, begin_idx_));
        auto end_it = make_iter<BulkPoint>( BulkPoint(elm_cache_map, element_patch_idx, end_idx_));
        return Range<BulkPoint>(bgn_it, end_it);
    }

protected:
    /// Index of data block according to subset in EvalPoints object.
    unsigned int subset_index_;
    uint begin_idx_;
    uint end_idx_;

};

/**
 * New Integral accessor class, replace of BulkIntegral, will be merged with BulkIntegral
 *
 * IN DEVELOPMENT
 */
template<unsigned int qdim>
class BulkIntegralAcc : public BulkIntegral, public FactoryBase<qdim>, public std::enable_shared_from_this<BulkIntegralAcc<qdim>> {
public:
    typedef BulkPoint PointType;
    typedef unsigned int MeshItem;

    /// Default constructor
    BulkIntegralAcc() : BulkIntegral() {}

    /// Constructor of bulk integral
    BulkIntegralAcc(std::shared_ptr<EvalPoints> eval_points, Quadrature *quad, PatchFEValues<3> *pfev, unsigned int i_subset)
     : BulkIntegral(eval_points, qdim, i_subset), FactoryBase<qdim>(pfev, quad)
    {
        this->fe_ = pfev->fe_dim<qdim>();
    }

    /// Destructor
    ~BulkIntegralAcc()
    {}

    /**
     * @brief Register the product of Jacobian determinant and the quadrature
     * weight at bulk quadrature points.
     *
     * @param quad Quadrature.
     */
    inline FeQ<Scalar> JxW()
    {
        return FeQ<Scalar>(this->template make_patch_op< Op::JxW<qdim, Op::BulkDomain, 3> >());
    }

	/// Create bulk accessor of coords entity
    inline FeQ<Vector> coords()
    {
        return FeQ<Vector>(this->template make_patch_op< Op::PtCoords<qdim, Op::BulkDomain, 3> >());
    }

//    inline ElQ<Tensor> jacobian(std::initializer_list<Quadrature *> quad_list)
//    {}

    /// Create bulk accessor of jac determinant entity
    inline ElQ<Scalar> determinant()
    {
        return ElQ<Scalar>( this->template make_patch_op< Op::JacDet<qdim, Op::BulkDomain, Op::BulkDomain, 3> >() );
    }

    /**
     * @brief Return the value of the @p function_no-th shape function at
     * the @p p bulk quadrature point.
     *
     * @param component_idx Number of the shape function.
     */
    inline FeQArray<Scalar> scalar_shape(uint component_idx = 0)
    {
        return this->template make_qarray<Scalar, Op::ScalarShape, Op::BulkDomain>(component_idx);
    }

    /**
     * @brief Return the value of the @p function_no-th vector shape function at
     * the @p p bulk quadrature point.
     *
     * @param component_idx Number of the shape function.
     */
    inline FeQArray<Vector> vector_shape(uint component_idx = 0)
    {
        return this->template make_qarray<Vector, Op::DispatchVectorShape, Op::BulkDomain>(component_idx);
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
        return this->template make_qarray<Vector, Op::GradScalarShape, Op::BulkDomain>(component_idx);
    }

    /**
     * @brief Return the value of the @p function_no-th gradient vector shape function
     * at the @p p bulk quadrature point.
     *
     * @param component_idx Number of the shape function.
     */
    inline FeQArray<Tensor> grad_vector_shape(uint component_idx=0)
    {
        return this->template make_qarray<Tensor, Op::DispatchGradVectorShape, Op::BulkDomain>(component_idx);
    }

    /**
     * @brief Return the value of the @p function_no-th vector symmetric gradient
     * at the @p p bulk quadrature point.
     *
     * @param component_idx Number of the shape function.
     */
    inline FeQArray<Tensor> vector_sym_grad(uint component_idx=0)
    {
        return this->template make_qarray<Tensor, Op::VectorSymGrad, Op::BulkDomain>(component_idx);
    }

    /**
     * @brief Return the value of the @p function_no-th vector divergence at
     * the @p p bulk quadrature point.
     *
     * @param component_idx Number of the shape function.
     */
    inline FeQArray<Scalar> vector_divergence(uint component_idx=0)
    {
        return this->template make_qarray<Scalar, Op::VectorDivergence, Op::BulkDomain>(component_idx);
    }
};

/**
 * Integral class of side points, allows assemblation of element - element fluxes.
 */
class EdgeIntegral : public BaseIntegral, public std::enable_shared_from_this<EdgeIntegral> {
public:
    /// Default constructor
	EdgeIntegral() : BaseIntegral()
    {
	    ASSERT_PERMANENT(false);
    }

    /// Constructor of edge integral
	EdgeIntegral(std::shared_ptr<EvalPoints> eval_points, unsigned int dim, uint i_subset);

    /// Destructor
    ~EdgeIntegral();

    /// Getter of n_sides
    inline unsigned int n_sides() const {
        return n_sides_;
    }

    /// Return index of data block according to subset in EvalPoints object
    inline int get_subset_idx() const {
        return subset_index_;
    }

    inline uint side_begin(const DHCellSide &cell_side) const {
        return begin_idx_ + cell_side.side_idx() * n_points_per_side_;
    }

    /// Returns range of side local points for appropriate cell side accessor
    inline Range< EdgePoint > points(const DHCellSide &cell_side, const ElementCacheMap *elm_cache_map) const {
        ASSERT_EQ(cell_side.dim(), dim_);
        //DebugOut() << "points per side: " << n_points_per_side_;
        uint element_patch_idx = elm_cache_map->position_in_cache(cell_side.element().idx());
        uint begin_idx = side_begin(cell_side);
        auto bgn_it = make_iter<EdgePoint>( EdgePoint(
                BulkPoint(elm_cache_map, element_patch_idx, 0), this, begin_idx));
        auto end_it = make_iter<EdgePoint>( EdgePoint(
                BulkPoint(elm_cache_map, element_patch_idx, n_points_per_side_), this, begin_idx));
        return Range<EdgePoint>(bgn_it, end_it);
    }


protected:
    unsigned int subset_index_;
    uint begin_idx_;

    /// Number of sides (value 0 indicates bulk set)
    unsigned int n_sides_;
    /// Number of points. TODO: pass this to the constructor, avoid extraction from the eval_points
    uint n_points_per_side_;

    friend class EvalPoints;
    friend class EdgePoint;
    friend class CouplingPoint;
    friend class BoundaryPoint;
    friend class CouplingIntegral;
    friend class BoundaryIntegral;
};

/**
 * New Integral accessor class, replace of EdgeIntegral, will be merged with EdgeIntegral
 *
 * IN DEVELOPMENT
 */
template<unsigned int qdim>
class EdgeIntegralAcc : public EdgeIntegral, public FactoryBase<qdim>, public std::enable_shared_from_this<BulkIntegralAcc<qdim>> {
public:
    typedef EdgePoint PointType;
    typedef DHCellSide MeshItem;

    /// Default constructor
    EdgeIntegralAcc() : EdgeIntegral() {}

    /// Constructor of bulk integral
    EdgeIntegralAcc(std::shared_ptr<EvalPoints> eval_points, Quadrature *quad, PatchFEValues<3> *pfev, unsigned int i_subset)
     : EdgeIntegral(eval_points, qdim, i_subset), FactoryBase<qdim>(pfev, quad)
    {
        this->fe_ = pfev->fe_dim<qdim>();
    }

    /// Destructor
    ~EdgeIntegralAcc()
    {}

    /// Same as BulkValues::JxW but register at side quadrature points.
    inline FeQ<Scalar> JxW()
    {
        return FeQ<Scalar>(this->template make_patch_op< Op::JxW<qdim, Op::SideDomain, 3> >());
    }

    /**
     * @brief Register the normal vector to a side at side quadrature points.
     *
     * @param quad Quadrature.
     */
	inline ElQ<Vector> normal_vector()
	{
        return ElQ<Vector>(this->template make_patch_op< Op::NormalVec<qdim, 3> >());
	}

	/// Create side accessor of coords entity
    inline FeQ<Vector> coords()
    {
        return FeQ<Vector>(this->template make_patch_op< Op::PtCoords<qdim, Op::SideDomain, 3> >());
    }

    /// Create bulk accessor of jac determinant entity
    inline ElQ<Scalar> determinant()
    {
        return ElQ<Scalar>(this->template make_patch_op< Op::JacDet<qdim, Op::SideDomain, Op::SideDomain, 3> >());
    }

    /// Same as BulkValues::scalar_shape but register at side quadrature points.
    inline FeQArray<Scalar> scalar_shape(uint component_idx = 0)
    {
        return this->template make_qarray<Scalar, Op::ScalarShape, Op::SideDomain>(component_idx);
    }

    /// Same as BulkValues::vector_shape but register at side quadrature points.
    inline FeQArray<Vector> vector_shape(uint component_idx = 0)
    {
        return this->template make_qarray<Vector, Op::DispatchVectorShape, Op::SideDomain>(component_idx);
    }

    /// Same as BulkValues::grad_scalar_shape but register at side quadrature points.
    inline FeQArray<Vector> grad_scalar_shape(uint component_idx=0)
    {
        return this->template make_qarray<Vector, Op::GradScalarShape, Op::SideDomain>(component_idx);
    }

    /**
     * @brief Return the value of the @p function_no-th gradient vector shape function
     * at the @p p bulk quadrature point.
     *
     * @param component_idx Number of the shape function.
     */
    inline FeQArray<Tensor> grad_vector_shape(uint component_idx=0)
    {
        return this->template make_qarray<Tensor, Op::DispatchGradVectorShape, Op::SideDomain>(component_idx);
    }

    /**
     * @brief Return the value of the @p function_no-th vector symmetric gradient
     * at the @p p side quadrature point.
     *
     * @param component_idx Number of the shape function.
     */
    inline FeQArray<Tensor> vector_sym_grad(uint component_idx=0)
    {
        return this->template make_qarray<Tensor, Op::VectorSymGrad, Op::SideDomain>(component_idx);
    }

    /**
     * @brief Return the value of the @p function_no-th vector divergence at
     * the @p p side quadrature point.
     *
     * @param component_idx Number of the shape function.
     */
    inline FeQArray<Scalar> vector_divergence(uint component_idx=0)
    {
        return this->template make_qarray<Scalar, Op::VectorDivergence, Op::SideDomain>(component_idx);
    }
};


/**
 * Integral class of neighbour points, allows assemblation of element - side fluxes.
 *
 * Dimension corresponds with element of higher dim.
 */
class CouplingIntegral : public BaseIntegral, public std::enable_shared_from_this<CouplingIntegral> {
public:
    /// Default constructor
	CouplingIntegral() : BaseIntegral() {}

    /// Constructor of ngh integral
	CouplingIntegral(std::shared_ptr<EdgeIntegral> edge_integral, std::shared_ptr<BulkIntegral> bulk_integral);

    /// Destructor
    ~CouplingIntegral();

    /// Return index of data block according to subset of higher dim in EvalPoints object
    inline int get_subset_high_idx() const {
        return edge_integral_->get_subset_idx();
    }

    /// Return index of data block according to subset of lower dim in EvalPoints object
    inline int get_subset_low_idx() const {
        return bulk_integral_->get_subset_idx();
    }

    inline uint bulk_begin() const {
        return eval_points_->subset_begin(dim_-1, bulk_integral_->get_subset_idx());
    }

    /// Returns range of side local points for appropriate cell side accessor
    inline Range< CouplingPoint > points(const DHCellSide &cell_side, const ElementCacheMap *elm_cache_map) const {
        ASSERT_EQ(cell_side.dim(), dim_);
        uint element_patch_idx = elm_cache_map->position_in_cache(cell_side.element().idx());
        uint begin_idx = edge_integral_->side_begin(cell_side);
        auto bgn_it = make_iter<CouplingPoint>( CouplingPoint(
                BulkPoint(elm_cache_map, element_patch_idx, 0), this, begin_idx) );
        auto end_it = make_iter<CouplingPoint>( CouplingPoint(
                BulkPoint(elm_cache_map, element_patch_idx, edge_integral_->n_points_per_side_), this, begin_idx) );;
        return Range<CouplingPoint>(bgn_it, end_it);
    }

private:
    /// Integral according to side subset part (element of higher dim) in EvalPoints object.
    std::shared_ptr<EdgeIntegral> edge_integral_;
    /// Integral according to bulk subset part (element of lower dim) in EvalPoints object.
    std::shared_ptr<BulkIntegral> bulk_integral_;

    friend class CouplingPoint;
};

template<unsigned int qdim>
class CouplingIntegralAcc : public CouplingIntegral, public FactoryBase<qdim>, public std::enable_shared_from_this<BoundaryIntegralAcc<1>> {
public:
    typedef CouplingPoint PointType;
    typedef DHCellSide MeshItem;

    /// Default constructor
    CouplingIntegralAcc() : CouplingIntegral() {}

    /// Constructor of bulk integral
    CouplingIntegralAcc(std::shared_ptr<EdgeIntegralAcc<qdim>> edge_integral, std::shared_ptr<BulkIntegralAcc<qdim-1>> bulk_integral,
            Quadrature *quad, PatchFEValues<3> *pfev)
    : CouplingIntegral(edge_integral, bulk_integral),
	  FactoryBase<qdim>(pfev, quad),
	  edge_integral_acc_(edge_integral), bulk_integral_acc_(bulk_integral)
    {
        this->fe_ = pfev->fe_dim<qdim>();
        fe_low_ = pfev->fe_dim<qdim-1>();
    }

    /// Destructor
    ~CouplingIntegralAcc()
    {
    	edge_integral_acc_.reset();
    	bulk_integral_acc_.reset();
    }

    /// Factory method. Same as previous but creates FE operation.
    template<class ValueType, template<unsigned int, class, unsigned int> class OpType>
    FeQJoin<ValueType> make_qjoin(uint component_idx = 0) {
        // element of lower dim (bulk points)
        auto fe_component_low = this->patch_fe_values_->fe_comp(fe_low_, component_idx);
        auto *low_dim_op = this->patch_fe_values_->template get< OpType<qdim-1, Op::BulkDomain, 3>, qdim-1 >(fe_component_low);
        auto *low_dim_zero_op = this->patch_fe_values_->template get< Op::OpZero<qdim-1, Op::BulkDomain, 3>, qdim-1 >(fe_component_low);

    	// element of higher dim (side points)
        auto fe_component_high = this->patch_fe_values_->fe_comp(this->fe_, component_idx);
        auto *high_dim_op = this->patch_fe_values_->template get< OpType<qdim, Op::SideDomain, 3>, qdim >(fe_component_high);
        auto *high_dim_zero_op = this->patch_fe_values_->template get< Op::OpZero<qdim, Op::SideDomain, 3>, qdim >(fe_component_high);

        ASSERT_EQ(fe_component_high->fe_type(), fe_component_low->fe_type()).error("Type of FiniteElement of low and high element must be same!\n");
        return FeQJoin<ValueType>(low_dim_op, high_dim_op, low_dim_zero_op, high_dim_zero_op);
    }

    /// Same as BulkValues::JxW but register at side quadrature points.
    inline FeQ<Scalar> JxW()
    {
        return edge_integral_acc_->JxW();
    }

    /**
     * @brief Register the normal vector to a side at side quadrature points.
     *
     * @param quad Quadrature.
     */
	inline ElQ<Vector> normal_vector()
	{
        return edge_integral_acc_->normal_vector();
	}

    inline FeQArray<Vector> vector_shape(uint component_idx = 0)
    {
        return edge_integral_acc_->vector_shape(component_idx);
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
    /// Integral according to higher dim (bulk) element subset part in EvalPoints object.
    std::shared_ptr<EdgeIntegralAcc<qdim>> edge_integral_acc_;
    /// Integral according to kower dim (boundary) element subset part in EvalPoints object.
    std::shared_ptr<BulkIntegralAcc<qdim-1>> bulk_integral_acc_;
    /// Holds FiniteEementt object of lower dimension
    std::shared_ptr< FiniteElement<qdim-1> > fe_low_;
};

/// Template specialization of previous class
template<>
class CouplingIntegralAcc<1> : public CouplingIntegral, public FactoryBase<1>, public std::enable_shared_from_this<CouplingIntegralAcc<1>> {
public:
    typedef CouplingPoint PointType;
    typedef DHCellSide MeshItem;

    /// Default constructor
    CouplingIntegralAcc() : CouplingIntegral() {}

    /// Constructor of bulk integral
    CouplingIntegralAcc(Quadrature *quad, PatchFEValues<3> *pfev)
    : CouplingIntegral(),
      FactoryBase<1>(pfev, quad)
    {
        this->fe_ = pfev->fe_dim<1>();
    }

    /// Destructor
    ~CouplingIntegralAcc()
    {}

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
};


/**
 * Integral class of boundary points, allows assemblation of fluxes between sides and neighbouring boundary elements.
 */
class BoundaryIntegral : public BaseIntegral, public std::enable_shared_from_this<BoundaryIntegral> {
public:
    /// Default constructor
    BoundaryIntegral() : BaseIntegral() {}

    /// Constructor of bulk subset
    BoundaryIntegral(std::shared_ptr<EdgeIntegral> edge_integral, std::shared_ptr<BulkIntegral> bulk_integral);

    /// Destructor
    ~BoundaryIntegral();

    /// Return index of data block according to subset of higher dim in EvalPoints object
    inline int get_subset_high_idx() const {
        return edge_integral_->get_subset_idx();
    }

    /// Return index of data block according to subset of lower dim (boundary) in EvalPoints object
    inline int get_subset_low_idx() const {
        return bulk_integral_->get_subset_idx();
    }

    inline uint bulk_begin() const {
      //  DebugOut().fmt("edge_begin: {} bdr_begin: {}",
      //          eval_points_->subset_begin(dim_, edge_integral_->get_subset_idx()),
      //          eval_points_->subset_begin(dim_-1, bulk_integral_->get_subset_idx()));
        return eval_points_->subset_begin(dim_-1, bulk_integral_->get_subset_idx());
    }

    /// Returns range of bulk local points for appropriate cell accessor
    inline Range< BoundaryPoint > points(const DHCellSide &cell_side, const ElementCacheMap *elm_cache_map) const {
        ASSERT_EQ(cell_side.dim(), dim_);
        uint element_patch_idx = elm_cache_map->position_in_cache(cell_side.element().idx());
        uint begin_idx = edge_integral_->side_begin(cell_side);
        auto bgn_it = make_iter<BoundaryPoint>( BoundaryPoint(
                BulkPoint(elm_cache_map, element_patch_idx, 0), this, begin_idx) );
        auto end_it = make_iter<BoundaryPoint>( BoundaryPoint(
                BulkPoint(elm_cache_map, element_patch_idx, edge_integral_->n_points_per_side_), this, begin_idx) );;
        return Range<BoundaryPoint>(bgn_it, end_it);
    }

protected:
    /// Integral according to higher dim (bulk) element subset part in EvalPoints object.
    std::shared_ptr<EdgeIntegral> edge_integral_;
    /// Integral according to kower dim (boundary) element subset part in EvalPoints object.
    std::shared_ptr<BulkIntegral> bulk_integral_;

    friend class BoundaryPoint;
};


template<unsigned int qdim>
class BoundaryIntegralAcc : public BoundaryIntegral, public FactoryBase<qdim>, public std::enable_shared_from_this<BoundaryIntegralAcc<qdim>> {
public:
    typedef BoundaryPoint PointType;
    typedef DHCellSide MeshItem;

    /// Default constructor
    BoundaryIntegralAcc() : BoundaryIntegral() {}

    /// Constructor of bulk integral
    BoundaryIntegralAcc(std::shared_ptr<EdgeIntegralAcc<qdim>> edge_integral, std::shared_ptr<BulkIntegralAcc<qdim-1>> bulk_integral,
            Quadrature *quad, PatchFEValues<3> *pfev)
    : BoundaryIntegral(edge_integral, bulk_integral),
	  FactoryBase<qdim>(pfev, quad),
	  edge_integral_acc_(edge_integral), bulk_integral_acc_(bulk_integral)
    {
        this->fe_ = pfev->fe_dim<qdim>();
    }

    /// Destructor
    ~BoundaryIntegralAcc()
    {
    	edge_integral_acc_.reset();
    	bulk_integral_acc_.reset();
    }

    /// Same as BulkValues::JxW but register at side quadrature points.
    inline FeQ<Scalar> JxW()
    {
        return edge_integral_acc_->JxW();
    }

    /**
     * @brief Register the normal vector to a side at side quadrature points.
     *
     * @param quad Quadrature.
     */
	inline ElQ<Vector> normal_vector()
	{
        return edge_integral_acc_->normal_vector();
	}

	/// Create side accessor of coords entity
    inline FeQ<Vector> coords()
    {
        return edge_integral_acc_->coords();
    }

    /// Create bulk accessor of jac determinant entity
    inline ElQ<Scalar> determinant()
    {
        return edge_integral_acc_->determinant();
    }

    /// Same as BulkValues::scalar_shape but register at side quadrature points.
    inline FeQArray<Scalar> scalar_shape(uint component_idx = 0)
    {
        return edge_integral_acc_->scalar_shape(component_idx);
    }

    /// Same as BulkValues::vector_shape but register at side quadrature points.
    inline FeQArray<Vector> vector_shape(uint component_idx = 0)
    {
        return edge_integral_acc_->vector_shape(component_idx);
    }

    /// Same as BulkValues::grad_scalar_shape but register at side quadrature points.
    inline FeQArray<Vector> grad_scalar_shape(uint component_idx=0)
    {
        return edge_integral_acc_->grad_scalar_shape(component_idx);
    }

    /**
     * @brief Return the value of the @p function_no-th gradient vector shape function
     * at the @p p bulk quadrature point.
     *
     * @param component_idx Number of the shape function.
     */
    inline FeQArray<Tensor> grad_vector_shape(uint component_idx=0)
    {
        return edge_integral_acc_->grad_vector_shape(component_idx);
    }

    /**
     * @brief Return the value of the @p function_no-th vector symmetric gradient
     * at the @p p side quadrature point.
     *
     * @param component_idx Number of the shape function.
     */
    inline FeQArray<Tensor> vector_sym_grad(uint component_idx=0)
    {
        return edge_integral_acc_->vector_sym_grad(component_idx);
    }

    /**
     * @brief Return the value of the @p function_no-th vector divergence at
     * the @p p side quadrature point.
     *
     * @param component_idx Number of the shape function.
     */
    inline FeQArray<Scalar> vector_divergence(uint component_idx=0)
    {
        return edge_integral_acc_->vector_divergence(component_idx);
    }

private:
    /// Integral according to higher dim (bulk) element subset part in EvalPoints object.
    std::shared_ptr<EdgeIntegralAcc<qdim>> edge_integral_acc_;
    /// Integral according to kower dim (boundary) element subset part in EvalPoints object.
    std::shared_ptr<BulkIntegralAcc<qdim-1>> bulk_integral_acc_;
};


#endif /* EVAL_SUBSET_HH_ */
