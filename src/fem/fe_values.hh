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
 * @file    fe_values.hh
 * @brief   Class FEValues calculates finite element data on the actual
 *          cells such as shape function values, gradients, Jacobian of
 *          the mapping from the reference cell etc.
 * @author  Jan Stebel
 */

#ifndef FE_VALUES_HH_
#define FE_VALUES_HH_

#include <string.h>                           // for memcpy
#include <algorithm>                          // for swap
#include <new>                                // for operator new[]
#include <string>                             // for operator<<
#include <vector>                             // for vector
#include "fem/element_values.hh"              // for ElementValues
#include "fem/fe_values_views.hh"             // for FEValuesViews
#include "mesh/ref_element.hh"                // for RefElement
#include "mesh/accessors.hh"
#include "fem/update_flags.hh"                // for UpdateFlags
#include "tools/mixed.hh"
#include "quadrature/quadrature_lib.hh"

class Quadrature;
template<unsigned int dim> class FiniteElement;

template<unsigned int dim> class MapScalar;
template<unsigned int dim> class MapPiola;
template<unsigned int dim> class MapContravariant;
template<unsigned int dim> class MapVector;
template<unsigned int dim> class MapTensor;
template<unsigned int dim> class MapSystem;

template<unsigned int spcedim> class FEValues;
template<unsigned int spcedim> class PatchFEValues;






template<class FV, unsigned int spacedim = 3>
class FEValuesBase
{
protected:
    // internal structure that stores all possible views
    // for scalar and vector-valued components of the FE
    struct ViewsCache {
        vector<FEValuesViews::Scalar<FV, spacedim> > scalars;
        vector<FEValuesViews::Vector<FV, spacedim> > vectors;
        vector<FEValuesViews::Tensor<FV, spacedim> > tensors;

        template<unsigned int DIM>
        void initialize(const FV &fv, const FiniteElement<DIM> &fe);
    };


public:
    /// Default constructor with postponed initialization.
    FEValuesBase();

    /**
	 * @brief Initialize structures and calculates cell-independent data.
	 *
	 * @param _quadrature The quadrature rule for the cell associated
     *                    to given finite element or for the cell side.
	 * @param _fe The finite element.
	 * @param _flags The update flags.
	 */
    template<unsigned int DIM>
    void initialize(Quadrature &_quadrature,
                    FiniteElement<DIM> &_fe,
                    UpdateFlags _flags);

    /**
     * @brief Allocates space for computed data.
     *
     * @param n_points    Number of quadrature points.
     * @param _fe         The finite element.
     * @param flags       The update flags.
     */
    template<unsigned int DIM>
    void allocate(unsigned int n_points,
                  FiniteElement<DIM> &_fe,
                  UpdateFlags flags);

    /**
     * @brief Returns the number of quadrature points.
     */
    inline unsigned int n_points() const
    { return n_points_; }

    /**
     * @brief Returns the number of shape functions.
     */
    inline unsigned int n_dofs() const
    {
        return n_dofs_;
    }

    /// Return dimension of reference space.
    inline unsigned int dim() const
    { return dim_; }

	/**
     * @brief Accessor to scalar values of multicomponent FE.
     * @param i Index of scalar component.
     */
	const FEValuesViews::Scalar<FV, spacedim> &scalar_view(unsigned int i) const
	{
      ASSERT_LT(i, views_cache_.scalars.size());
      return views_cache_.scalars[i];
    }

    /**
     * @brief Accessor to vector values of multicomponent FE.
     * @param i Index of first vector component.
     */
    const FEValuesViews::Vector<FV, spacedim> &vector_view(unsigned int i) const
    {
      ASSERT_LT(i, views_cache_.vectors.size());
      return views_cache_.vectors[i];
    }

    /**
     * @brief Accessor to tensor values of multicomponent FE.
     * @param i Index of first tensor component.
     */
    const FEValuesViews::Tensor<FV, spacedim> &tensor_view(unsigned int i) const
    {
      ASSERT_LT(i, views_cache_.tensors.size());
      return views_cache_.tensors[i];
    }

protected:
    /// Structure for storing the precomputed finite element data.
    class FEInternalData
    {
    public:

        FEInternalData(unsigned int np, unsigned int nd);

        /// Create a new instance of FEInternalData for a FESystem component or subvector.
        FEInternalData(const FEInternalData &fe_system_data,
                       const std::vector<unsigned int> &dof_indices,
                       unsigned int first_component_idx,
                       unsigned int ncomponents = 1);

        /**
         * @brief Precomputed values of basis functions at the quadrature points.
         *
         * Dimensions:   (no. of quadrature points)
         *             x (no. of dofs)
         *             x (no. of components in ref. cell)
         */
        std::vector<std::vector<arma::vec> > ref_shape_values;

        /**
         * @brief Precomputed gradients of basis functions at the quadrature points.
         *
         * Dimensions:   (no. of quadrature points)
         *             x (no. of dofs)
         *             x ((dim_ of. ref. cell)x(no. of components in ref. cell))
         */
        std::vector<std::vector<arma::mat> > ref_shape_grads;

        /// Number of quadrature points.
        unsigned int n_points;

        /// Number of dofs (shape functions).
        unsigned int n_dofs;
    };


    /// Initialize vectors declared separately in descendants.
    virtual void allocate_in() =0;

    /// Initialize ElementValues separately in descendants.
    virtual void initialize_in(Quadrature &, unsigned int) =0;

    /// Precompute finite element data on reference element.
    template<unsigned int DIM>
    std::shared_ptr<typename FEValuesBase<FV, spacedim>::FEInternalData> init_fe_data(const FiniteElement<DIM> &fe, const Quadrature &q);

    /// Dimension of reference space.
    unsigned int dim_;

    /// Number of integration points.
    unsigned int n_points_;

    /// Number of finite element dofs.
    unsigned int n_dofs_;

    /// Type of finite element (scalar, vector, tensor).
    FEType fe_type_;

    /// Dof indices of FESystem sub-elements.
    std::vector<std::vector<unsigned int>> fe_sys_dofs_;

    /// Numbers of components of FESystem sub-elements in reference space.
    std::vector<unsigned int> fe_sys_n_components_;

    /// Numbers of components of FESystem sub-elements in real space.
    std::vector<unsigned int> fe_sys_n_space_components_;

    /// Flags that indicate which finite element quantities are to be computed.
    UpdateFlags update_flags;

    /// Vector of FEValues for sub-elements of FESystem.
    std::vector<FEValues<spacedim>> fe_values_vec;

    /// Number of components of the FE.
    unsigned int n_components_;

    /// Auxiliary storage of FEValuesViews accessors.
    ViewsCache views_cache_;

    /// Precomputed finite element data.
    std::shared_ptr<FEInternalData> fe_data_;

    /// Precomputed FE data (shape functions on reference element) for all side quadrature points.
    std::vector<shared_ptr<FEInternalData> > side_fe_data_;

    /// Helper object, we need its for ViewsCache initialization.
    const FV *fv_;

};





/**
 * @brief Calculates finite element data on the actual cell.
 *
 * FEValues takes care of the calculation of finite element data on
 * the actual cell or on its side, (values of shape functions at quadrature
 * points, gradients of shape functions, Jacobians of the mapping
 * from the reference cell etc.).
 * @param spacedim Dimension of the Euclidean space where the actual
 *                 cell lives.
 */
template<unsigned int spacedim = 3>
class FEValues : public FEValuesBase<FEValues<spacedim>, spacedim>
{
public:

    /// Default constructor with postponed initialization.
    FEValues();

    /// Constructor with initialization of data structures
    /// (see initialize() for description of parameters).
    template<unsigned int DIM>
    FEValues(Quadrature &_quadrature,
             FiniteElement<DIM> &_fe,
             UpdateFlags _flags)
    : FEValues()
    {
        this->initialize(_quadrature, _fe, _flags);
    }


    /// Correct deallocation of objects created by 'initialize' methods.
    ~FEValues();


    /**
     * @brief Update cell-dependent data (gradients, Jacobians etc.)
     *
     * @param cell The actual cell.
     */
    void reinit(const ElementAccessor<spacedim> &cell);

    /**
	 * @brief Update cell side-dependent FE data (values, gradients).
	 *
	 * @param cell_side Accessor to cell side.
	 */
    void reinit(const Side &cell_side);
    
    /**
     * @brief Return the value of the @p function_no-th shape function at
     * the @p point_no-th quadrature point.
     *
     * @param function_no Number of the shape function.
     * @param point_no Number of the quadrature point.
     */
    inline double shape_value(const unsigned int function_no, const unsigned int point_no) const
    {
        ASSERT_LT(function_no, this->n_dofs_);
        ASSERT_LT(point_no, this->n_points_);
        return shape_values_[point_no][function_no];
    }


    /**
     * @brief Return the gradient of the @p function_no-th shape function at
     * the @p point_no-th quadrature point.
     *
     * @param function_no Number of the shape function.
     * @param point_no Number of the quadrature point.
     */
    inline arma::vec::fixed<spacedim> shape_grad(const unsigned int function_no, const unsigned int point_no) const
	{
        ASSERT_LT(function_no, this->n_dofs_);
        ASSERT_LT(point_no, this->n_points_);
        return shape_gradients_[point_no][function_no];
    }

    /**
     * @brief Return the value of the @p function_no-th shape function at
     * the @p point_no-th quadrature point.
     *
     * For vectorial finite elements.
     *
     * @param function_no Number of the shape function.
     * @param point_no Number of the quadrature point.
     */
    inline double shape_value_component(const unsigned int function_no,
                                        const unsigned int point_no, 
                                        const unsigned int comp) const
    {
        ASSERT_LT(function_no, this->n_dofs_);
        ASSERT_LT(point_no, this->n_points_);
        ASSERT_LT(comp, this->n_components_);
        return shape_values_[point_no][function_no*this->n_components_+comp];
    }

    /**
     * @brief Return the gradient of the @p function_no-th shape function at
     * the @p point_no-th quadrature point.
     *
     * For vectorial finite elements.
     *
     * @param function_no Number of the shape function.
     * @param point_no Number of the quadrature point.
     */
    arma::vec::fixed<spacedim> shape_grad_component(const unsigned int function_no,
                                                           const unsigned int point_no,
                                                           const unsigned int comp) const;

    /**
     * @brief Return the relative volume change of the cell (Jacobian determinant).
     *
     * If dim_==spacedim then the sign may be negative, otherwise the
     * result is a positive number.
     *
     * @param point_no Number of the quadrature point.
     */
    inline double determinant(const unsigned int point_no)
    {
        ASSERT_LT(point_no, this->n_points_);
        return elm_values_->determinant(point_no);
    }

    /**
     * @brief Return the product of Jacobian determinant and the quadrature
     * weight at given quadrature point.
     *
     * @param point_no Number of the quadrature point.
     */
    inline double JxW(const unsigned int point_no)
    {
        ASSERT_LT(point_no, this->n_points_);
        // TODO: This is temporary solution to distinguish JxW on element and side_JxW on side.
        // In future we should call the appropriate method in elm_values_.
        return (elm_values_->cell().is_valid()) ? elm_values_->JxW(point_no) : elm_values_->side_JxW(point_no);
    }

    /**
     * @brief Return coordinates of the quadrature point in the actual cell system.
     *
     * @param point_no Number of the quadrature point.
     */
    inline arma::vec::fixed<spacedim> point(const unsigned int point_no)
    {
        ASSERT_LT(point_no, this->n_points_);
        return elm_values_->point(point_no);
    }

    /**
	 * @brief Return coordinates of all quadrature points in the actual cell system.
	 *
	 */
	inline const Armor::array &point_list() const
	{
	    return elm_values_->point_list();
	}


    /**
     * @brief Returns the normal vector to a side at given quadrature point.
     *
     * @param point_no Number of the quadrature point.
     */
	inline arma::vec::fixed<spacedim> normal_vector(unsigned int point_no)
	{
        ASSERT_LT(point_no, this->n_points_);
	    return elm_values_->normal_vector(point_no);
	}
	

protected:
    
    /**
     * @brief Computes the shape function values and gradients on the actual cell
     * and fills the FEValues structure.
     *
     * @param fe_data Precomputed finite element data.
     */
    void fill_data(const ElementValues<spacedim> &elm_values, const typename FEValuesBase<FEValues<spacedim>, spacedim>::FEInternalData &fe_data);
    
    /**
     * @brief Computes the shape function values and gradients on the actual cell
     * and fills the FEValues structure. Specialized variant of previous method for
     * different FETypes given by template parameter.
     */
    template<class MapType>
    void fill_data_specialized(const ElementValues<spacedim> &elm_values, const typename FEValuesBase<FEValues<spacedim>, spacedim>::FEInternalData &fe_data);


    /// Implement @p FEValuesBase::allocate_in
    void allocate_in() override;

    /// Implement @p FEValuesBase::initialize_in
    void initialize_in(Quadrature &q, unsigned int dim) override;

    /// Shape functions evaluated at the quadrature points.
    std::vector<std::vector<double> > shape_values_;

    /// Gradients of shape functions evaluated at the quadrature points.
    /// Each row of the matrix contains the gradient of one shape function.
    std::vector<std::vector<arma::vec::fixed<spacedim> > > shape_gradients_;

    /// Auxiliary object for calculation of element-dependent data.
    std::shared_ptr<ElementValues<spacedim> > elm_values_;
    
    friend class MapScalar<spacedim>;
    friend class MapPiola<spacedim>;
    friend class MapContravariant<spacedim>;
    friend class MapVector<spacedim>;
    friend class MapTensor<spacedim>;
    friend class MapSystem<spacedim>;
};






template<unsigned int spacedim = 3>
class PatchFEValues : public FEValuesBase<PatchFEValues<spacedim>, spacedim> {
public:
    /// Constructor, set size of ElementData vector (maximal capacity)
	PatchFEValues(unsigned int max_size);

    /// Set new value of used_size_
    void resize(unsigned int new_size);

	inline unsigned int used_size() const {
	    return used_size_;
	}

	inline unsigned int max_size() const {
	    return element_data_.size();
	}

    /**
     * @brief Return the value of the @p function_no-th shape function at
     * the @p point_no-th quadrature point.
     *
     * For vectorial finite elements.
     *
     * @param function_no Number of the shape function.
     * @param point_no Number of the quadrature point.
     */
    inline double shape_value_component(const unsigned int function_no,
                                        const unsigned int point_no,
                                        const unsigned int comp) const
    {
        ASSERT_LT(function_no, this->n_dofs_);
        ASSERT_LT(point_no, this->n_points_);
        ASSERT_LT(comp, this->n_components_);
        return element_data_[patch_cell_idx_].shape_values_[point_no][function_no*this->n_components_+comp];
    }

    /**
     * @brief Return the gradient of the @p function_no-th shape function at
     * the @p point_no-th quadrature point.
     *
     * For vectorial finite elements.
     *
     * @param function_no Number of the shape function.
     * @param point_no Number of the quadrature point.
     */
    arma::vec::fixed<spacedim> shape_grad_component(const unsigned int function_no,
                                                           const unsigned int point_no,
                                                           const unsigned int comp) const;

protected:
    class ElementFEData
    {
    public:
        ElementFEData() {}

        /// Shape functions evaluated at the quadrature points.
        std::vector<std::vector<double> > shape_values_;

        /// Gradients of shape functions evaluated at the quadrature points.
        /// Each row of the matrix contains the gradient of one shape function.
        std::vector<std::vector<arma::vec::fixed<spacedim> > > shape_gradients_;

        /// Auxiliary object for calculation of element-dependent data.
        std::shared_ptr<ElementValues<spacedim> > elm_values_;

    };

    /// Implement @p FEValuesBase::allocate_in
    void allocate_in() override;

    /// Implement @p FEValuesBase::initialize_in
    void initialize_in(Quadrature &q, unsigned int dim) override;

    /// Patch index of processed element.
    unsigned int patch_cell_idx_;

    /// Index of patch of processed element.
    unsigned int side_idx_;

    /// Data of elements on patch
    std::vector<ElementFEData> element_data_;

    /// Number of elements on patch. Must be less or equal to size of element_data vector
    unsigned int used_size_;

};



std::vector<FEValues<3>> mixed_fe_values(
        QGauss::array &quadrature,
        MixedPtr<FiniteElement> fe,
        UpdateFlags flags);







#endif /* FE_VALUES_HH_ */
