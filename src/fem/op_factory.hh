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
 * @file    op_factory.hh
 * @brief   Declares top level factory classes of FE operations.
 * @author  David Flanderka
 */

#ifndef OP_FACTORY_HH_
#define OP_FACTORY_HH_

#include "fem/eigen_tools.hh"
#include "fem/op_function.hh"
#include "fem/op_accessors_impl.hh"

template<unsigned int spacedim> class PatchFEValues;


template<unsigned int dim>
class BaseValues
{
protected:
	// Default constructor
	BaseValues(PatchFEValues<3> &pfev) : patch_fe_values_(pfev)
	{}

    /// Return FiniteElement of \p component_idx for FESystem or \p fe for other types
    template<unsigned int FE_dim>
    std::shared_ptr<FiniteElement<FE_dim>> fe_comp(std::shared_ptr< FiniteElement<FE_dim> > fe, uint component_idx) {
        if (fe->fe_type() == FEMixedSystem) {
            FESystem<FE_dim> *fe_sys = dynamic_cast<FESystem<FE_dim>*>( fe.get() );
            return fe_sys->fe()[component_idx];
        } else {
            ASSERT_EQ(component_idx, 0).warning("Non-zero component_idx can only be used for FESystem.");
            return fe;
        }
    }

    /// Factory method. Creates operation of given OpType.
    template<class OpType>
    PatchOp<3> *make_patch_op() {
    	return patch_fe_values_.get< OpType, dim >();
    }

    /// Factory method. Same as previous but creates FE operation.
    template<class OpType>
    PatchOp<3> *make_patch_op(std::shared_ptr<FiniteElement<dim>> fe) {
    	return patch_fe_values_.get< OpType, dim >(fe);
    }

    /// Factory method. Same as previous but creates FE operation.
    template<class ValueType, template<unsigned int, class, unsigned int> class OpType, class Domain>
    FeQArray<ValueType> make_qarray(uint component_idx = 0) {
    	std::shared_ptr<FiniteElement<dim>> fe_component = this->fe_comp(fe_, component_idx);
    	return FeQArray<ValueType>(this->template make_patch_op< OpType<dim, Domain, 3> >(fe_component));
    }

    /// Factory method. Same as previous but creates FE operation.
    template<class ValueType, template<unsigned int, unsigned int> class OpType>
    FeQArray<ValueType> make_qarray2(uint component_idx = 0) {
    	std::shared_ptr<FiniteElement<dim>> fe_component = this->fe_comp(fe_, component_idx);
    	return FeQArray<ValueType>(this->template make_patch_op< OpType<dim, 3> >(fe_component));
    }

    PatchFEValues<3> &patch_fe_values_;
    std::shared_ptr< FiniteElement<dim> > fe_;
};

template<unsigned int dim>
class BulkValues : public BaseValues<dim>
{
public:
	/// Constructor
	BulkValues(PatchFEValues<3> &pfev, MixedPtr<FiniteElement> fe)
	: BaseValues<dim>(pfev) {
	    this->fe_ = fe[Dim<dim>{}];
	}

    /**
     * @brief Register the product of Jacobian determinant and the quadrature
     * weight at bulk quadrature points.
     *
     * @param quad Quadrature.
     */
    inline FeQ<Scalar> JxW()
    {
        return FeQ<Scalar>(this->template make_patch_op< Op::JxW<dim, Op::BulkDomain, 3> >());
    }

	/// Create bulk accessor of coords entity
    inline FeQ<Vector> coords()
    {
        return FeQ<Vector>(this->template make_patch_op< Op::PtCoords<dim, Op::BulkDomain, 3> >());
    }

//    inline ElQ<Tensor> jacobian(std::initializer_list<Quadrature *> quad_list)
//    {}

    /// Create bulk accessor of jac determinant entity
    inline ElQ<Scalar> determinant()
    {
        return ElQ<Scalar>( this->template make_patch_op< Op::JacDet<dim, Op::BulkDomain, Op::BulkDomain, 3> >() );
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


template<unsigned int dim>
class SideValues : public BaseValues<dim>
{
public:
	/// Constructor
	SideValues(PatchFEValues<3> &pfev, MixedPtr<FiniteElement> fe)
	: BaseValues<dim>(pfev) {
	    this->fe_ = fe[Dim<dim>{}];
	}

    /// Same as BulkValues::JxW but register at side quadrature points.
    inline FeQ<Scalar> JxW()
    {
        return FeQ<Scalar>(this->template make_patch_op< Op::JxW<dim, Op::SideDomain, 3> >());
    }

    /**
     * @brief Register the normal vector to a side at side quadrature points.
     *
     * @param quad Quadrature.
     */
	inline ElQ<Vector> normal_vector()
	{
        return ElQ<Vector>(this->template make_patch_op< Op::Side::Pt::OpNormalVec<dim, 3> >());
	}

	/// Create side accessor of coords entity
    inline FeQ<Vector> coords()
    {
        return FeQ<Vector>(this->template make_patch_op< Op::PtCoords<dim, Op::SideDomain, 3> >());
    }

    /// Create bulk accessor of jac determinant entity
    inline ElQ<Scalar> determinant()
    {
        return ElQ<Scalar>(this->template make_patch_op< Op::JacDet<dim, Op::SideDomain, Op::SideDomain, 3> >());
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


template<unsigned int dim>
class JoinValues : public BaseValues<dim>
{
public:
	/// Constructor
	JoinValues(PatchFEValues<3> &pfev, MixedPtr<FiniteElement> fe)
	: BaseValues<dim>(pfev) {
	    fe_high_dim_ = fe[Dim<dim>{}];
	    fe_low_dim_ = fe[Dim<dim-1>{}];
	}

    inline FeQJoin<Scalar> scalar_join_shape(uint component_idx = 0)
    {
    	// element of lower dim (bulk points)
        auto fe_component_low = this->fe_comp(fe_low_dim_, component_idx);
        ASSERT_EQ(fe_component_low->fe_type(), FEType::FEScalar).error("Type of FiniteElement of scalar_shape accessor must be FEScalar!\n");
        auto *low_dim_op = this->patch_fe_values_.template get< Op::ScalarShape<dim-1, Op::BulkDomain, 3>, dim-1 >(fe_component_low);
        auto *low_dim_zero_op = this->patch_fe_values_.template get< Op::OpZero<dim-1, Op::BulkDomain, 3>, dim-1 >(fe_component_low);

    	// element of higher dim (side points)
        auto fe_component_high = this->fe_comp(fe_high_dim_, component_idx);
        ASSERT_EQ(fe_component_high->fe_type(), FEType::FEScalar).error("Type of FiniteElement of scalar_shape accessor must be FEScalar!\n");
        auto *high_dim_op = this->template make_patch_op< Op::ScalarShape<dim, Op::SideDomain, 3> >(fe_component_high);
        auto *high_dim_zero_op = this->template make_patch_op< Op::OpZero<dim, Op::SideDomain, 3> >(fe_component_high);

        return FeQJoin<Scalar>(low_dim_op, high_dim_op, low_dim_zero_op, high_dim_zero_op);
    }

    inline FeQJoin<Vector> vector_join_shape(uint component_idx = 0)
    {
    	// element of lower dim (bulk points)
        auto fe_component_low = this->fe_comp(fe_low_dim_, component_idx);
        auto *low_dim_op = this->patch_fe_values_.template get< Op::DispatchVectorShape<dim-1, Op::BulkDomain, 3>, dim-1 >(fe_component_low);
        auto *low_dim_zero_op = this->patch_fe_values_.template get< Op::OpZero<dim-1, Op::BulkDomain, 3>, dim-1 >(fe_component_low);

        // element of higher dim (side points)
        auto fe_component_high = this->fe_comp(fe_high_dim_, component_idx);
        auto *high_dim_op = this->template make_patch_op< Op::DispatchVectorShape<dim, Op::SideDomain, 3> >(fe_component_high);
        auto *high_dim_zero_op = this->template make_patch_op< Op::OpZero<dim, Op::SideDomain, 3> >(fe_component_high);

        ASSERT_EQ(fe_component_high->fe_type(), fe_component_low->fe_type()).error("Type of FiniteElement of low and high element must be same!\n");
        return FeQJoin<Vector>(low_dim_op, high_dim_op, low_dim_zero_op, high_dim_zero_op);
    }

    inline FeQJoin<Tensor> gradient_vector_join_shape(uint component_idx = 0)
    {
    	// element of lower dim (bulk points)
        auto fe_component_low = this->fe_comp(fe_low_dim_, component_idx);
        auto *low_dim_op = this->patch_fe_values_.template get< Op::DispatchGradVectorShape<dim-1, Op::BulkDomain, 3>, dim-1 >(fe_component_low);
        auto *low_dim_zero_op = this->patch_fe_values_.template get< Op::OpZero<dim-1, Op::BulkDomain, 3>, dim-1 >(fe_component_low);

        // element of higher dim (side points)
        auto fe_component_high = this->fe_comp(fe_high_dim_, component_idx);
        auto *high_dim_op = this->template make_patch_op< Op::DispatchGradVectorShape<dim, Op::SideDomain, 3> >(fe_component_high);
        auto *high_dim_zero_op = this->template make_patch_op< Op::OpZero<dim, Op::SideDomain, 3> >(fe_component_high);

        ASSERT_EQ(fe_component_high->fe_type(), fe_component_low->fe_type()).error("Type of FiniteElement of low and high element must be same!\n");
        return FeQJoin<Tensor>(low_dim_op, high_dim_op, low_dim_zero_op, high_dim_zero_op);
    }

private:
    std::shared_ptr< FiniteElement<dim> > fe_high_dim_;
    std::shared_ptr< FiniteElement<dim-1> > fe_low_dim_;
};

/// Template specialization of dim = 1
template <>
class JoinValues<1> : public BaseValues<1>
{
public:
	/// Constructor
	JoinValues(PatchFEValues<3> &pfev, FMT_UNUSED MixedPtr<FiniteElement> fe)
	: BaseValues<1>(pfev) {}

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


#endif /* OP_FACTORY_HH_ */
