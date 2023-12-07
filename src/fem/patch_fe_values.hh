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
 * @file    patch_fe_values.hh
 * @brief   Class FEValues calculates finite element data on the actual
 *          cells such as shape function values, gradients, Jacobian of
 *          the mapping from the reference cell etc.
 * @author  Jan Stebel, David Flanderka
 */

#ifndef PATCH_FE_VALUES_HH_
#define PATCH_FE_VALUES_HH_


#include <string.h>                           // for memcpy
#include <algorithm>                          // for swap
#include <new>                                // for operator new[]
#include <string>                             // for operator<<
#include <vector>                             // for vector
#include "fem/element_values.hh"              // for ElementValues
#include "fem/fe_values.hh"                   // for FEValuesBase
#include "fem/fe_values_views.hh"             // for FEValuesViews
#include "mesh/ref_element.hh"                // for RefElement
#include "mesh/accessors.hh"
#include "fem/update_flags.hh"                // for UpdateFlags
#include "quadrature/quadrature_lib.hh"
#include "fields/eval_subset.hh"

template<unsigned int spacedim> class PatchFEValues;



using Scalar = double;
using Vector = arma::vec3;
using Tensor = arma::mat33;

template <class ValueType>
class ElQ {
public:
    /// Forbidden default constructor
    ElQ() = delete;

    /// Constructor
    ElQ(PatchFEValues<3> *fe_values, unsigned int begin)
    : fe_values_(fe_values), begin_(begin) {}

    ValueType operator()(FMT_UNUSED const BulkPoint &point);

private:
    // attributes:
    PatchFEValues<3> *fe_values_;
    unsigned int begin_;    /// Index of the first component of the Quantity. Size is given by ValueType
};


template <class ValueType>
class FeQ {
public:
    /// Forbidden default constructor
    FeQ() = delete;

    // Class similar to current FeView
    FeQ(PatchFEValues<3> *fe_values, unsigned int begin)
    : fe_values_(fe_values), begin_(begin) {}


    ValueType operator()(FMT_UNUSED unsigned int shape_idx, FMT_UNUSED const BulkPoint &point);

    // Implementation for EdgePoint, SidePoint, and JoinPoint shoud have a common implementation
    // resolving to side values

private:
    // attributes:
    PatchFEValues<3> *fe_values_;
    unsigned int begin_;    /// Index of the first component of the Quantity. Size is given by ValueType
};


template<unsigned int spacedim = 3>
class PatchFEValues {
public:

    PatchFEValues(unsigned int n_quad_points)
    : dim_fe_vals_({DimPatchFEValues(n_quad_points), DimPatchFEValues(n_quad_points), DimPatchFEValues(n_quad_points)}),
	  n_columns_(0) {}

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
                    UpdateFlags _flags)
    {
        dim_fe_vals_[DIM-1].initialize(_quadrature, _fe, _flags);
    }

    /// Reinit data.
    void reinit(std::array<PatchElementsList, 4> patch_elements)
    {
        for (unsigned int i=0; i<3; ++i) {
            dim_fe_vals_[i].reinit(patch_elements[i+1]);
        }
    }

    /**
     * @brief Return the product of Jacobian determinant and the quadrature
     * weight at given quadrature point.
     *
     * @param quad_list List of quadratures.
     */
    inline ElQ<Scalar> JxW(FMT_UNUSED std::vector<Quadrature *> quad_list)
    {
        uint begin = this->n_columns_;
        n_columns_++; // scalar needs one column
        // TODO store to map?? JxW will be pre-computed to column 'begin'.
        return ElQ<Scalar>(this, begin);
    }

    /**
     * @brief Returns the normal vector to a side at given quadrature point.
     *
     * @param quad_list List of quadratures.
     */
	inline ElQ<Vector> normal_vector(FMT_UNUSED std::vector<Quadrature *> quad_list)
	{
        uint begin = this->n_columns_;
        n_columns_ += 3; // Vector needs 3 columns
        // TODO store to map?? shape_value will be pre-computed to column 'begin'.
        return ElQ<Vector>(this, begin);
	}

    /**
     * @brief Return the value of the @p function_no-th shape function at
     * the @p p quadrature point.
     *
     * @param quad_list List of quadratures.
     * @param function_no Number of the shape function.
     */
    inline FeQ<Scalar> scalar_shape(FMT_UNUSED std::vector<Quadrature *> quad_list, unsigned int n_comp)
    {
        uint begin = this->n_columns_;
        n_columns_ += n_comp; // scalar needs one column x n_comp
        // TODO store to map?? shape_value will be pre-computed to column 'begin'.
        return FeQ<Scalar>(this, begin);
    }

    inline FeQ<Vector> grad_scalar_shape(FMT_UNUSED std::vector<Quadrature *> quad_list, unsigned int n_comp)
    {
        uint begin = this->n_columns_;
        n_columns_ += 3 * n_comp; // scalar needs one column x n_comp
        // TODO store to map?? shape_value will be pre-computed to column 'begin'.
        return FeQ<Vector>(this, begin);
    }

private:
    enum MeshObjectType {
        ElementFE = 0,
		SideFE = 1
    };


//    /// Structure for storing the precomputed finite element data.
//    class FEInternalData
//    {
//    public:
//
//        FEInternalData(unsigned int np, unsigned int nd);
//
//        /// Create a new instance of FEInternalData for a FESystem component or subvector.
//        FEInternalData(const FEInternalData &fe_system_data,
//                       const std::vector<unsigned int> &dof_indices,
//                       unsigned int first_component_idx,
//                       unsigned int ncomps = 1);
//
//        /**
//         * @brief Precomputed values of basis functions at the quadrature points.
//         *
//         * Dimensions:   (no. of quadrature points)
//         *             x (no. of dofs)
//         *             x (no. of components in ref. cell)
//         */
//        std::vector<std::vector<arma::vec> > ref_shape_values;
//
//        /**
//         * @brief Precomputed gradients of basis functions at the quadrature points.
//         *
//         * Dimensions:   (no. of quadrature points)
//         *             x (no. of dofs)
//         *             x ((dim_ of. ref. cell)x(no. of components in ref. cell))
//         */
//        std::vector<std::vector<arma::mat> > ref_shape_grads;
//
//        /// Number of quadrature points.
//        unsigned int n_points;
//
//        /// Number of dofs (shape functions).
//        unsigned int n_dofs;
//    };


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


    /// Subobject holds FE data of one dimension (0,1,2,3)
    class DimPatchFEValues {
    public:
        /// Constructor
        DimPatchFEValues(unsigned int max_size=0)
        : used_size_(0), max_n_elem_(max_size) {}


    	inline unsigned int used_size() const {
    	    return used_size_;
    	}

    	inline unsigned int max_size() const {
    	    return element_data_.size();
    	}

        void reinit(PatchElementsList patch_elements);

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
        void allocate(Quadrature &_quadrature,
                      FiniteElement<DIM> &_fe,
                      UpdateFlags flags);

        /// Precompute finite element data on reference element.
        template<unsigned int DIM>
        std::shared_ptr<FEInternalData> init_fe_data(const FiniteElement<DIM> &fe, const Quadrature &q);

        /**
         * @brief Computes the shape function values and gradients on the actual cell
         * and fills the FEValues structure.
         *
         * @param fe_data Precomputed finite element data.
         */
        void fill_data(const ElementValues<spacedim> &elm_values, const FEInternalData &fe_data);

        /**
         * @brief Computes the shape function values and gradients on the actual cell
         * and fills the FEValues structure. Specialized variant of previous method for
         * different FETypes given by template parameter.
         */
        template<class MapType>
        void fill_data_specialized(const ElementValues<spacedim> &elm_values, const FEInternalData &fe_data);

        /**
         * Temporary function. Use in fill_data.
         * Set shape value @p val of the @p i_point and @p i_func_comp.
         */
        inline void set_shape_value(unsigned int i_point, unsigned int i_func_comp, double val)
        {
            element_data_[patch_data_idx_].shape_values_[i_point][i_func_comp] = val;
        }

        /**
         * Temporary function. Use in fill_data.
         * Set shape gradient @p val of the @p i_point and @p i_func_comp.
         */
        inline void set_shape_gradient(unsigned int i_point, unsigned int i_func_comp, arma::vec::fixed<spacedim> val)
        {
            element_data_[patch_data_idx_].shape_gradients_[i_point][i_func_comp] = val;
        }

        /**
         * Temporary function. Use in fill_data.
         * @brief Return the value of the @p function_no-th shape function at
         * the @p point_no-th quadrature point.
         *
         * @param function_no Number of the shape function.
         * @param point_no Number of the quadrature point.
         */
        inline double shape_value(const unsigned int function_no, const unsigned int point_no) const
        {
            // TODO: obsolete method, will be replaced with following two methods.
            ASSERT_LT(function_no, this->n_dofs_);
            ASSERT_LT(point_no, this->n_points_);
            return element_data_[patch_data_idx_].shape_values_[point_no][function_no];
        }

        /**
         * Temporary function. Use in fill_data.
         * @brief Return the gradient of the @p function_no-th shape function at
         * the @p point_no-th quadrature point.
         *
         * @param function_no Number of the shape function.
         * @param point_no Number of the quadrature point.
         */
        inline arma::vec::fixed<spacedim> shape_grad(const unsigned int function_no, const unsigned int point_no) const
    	{
            // TODO: obsolete method, will be replaced with following two methods.
            ASSERT_LT(function_no, this->n_dofs_);
            ASSERT_LT(point_no, this->n_points_);
            return element_data_[patch_data_idx_].shape_gradients_[point_no][function_no];
        }

        /// Set size of ElementFEData. Important: Use only during the initialization of FESystem !
        void resize(unsigned int max_size) {
            element_data_.resize(max_size);
        }


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
        std::vector<PatchFEValues<spacedim>::DimPatchFEValues> fe_values_vec;

        /// Number of components of the FE.
        unsigned int n_components_;

//        /// Auxiliary storage of FEValuesViews accessors.
//        ViewsCache views_cache_;

        /// Precomputed finite element data.
        std::shared_ptr<FEInternalData> fe_data_;

        /// Precomputed FE data (shape functions on reference element) for all side quadrature points.
        std::vector<shared_ptr<FEInternalData> > side_fe_data_;

        /// Patch index of processed element / side.
        unsigned int patch_data_idx_;

        /// Map of element patch indexes to element_data_.
        std::map<unsigned int, unsigned int> element_patch_map_;

        /// Data of elements / sides on patch
        std::vector<ElementFEData> element_data_;

        /// Number of elements / sides on patch. Must be less or equal to size of element_data vector
        unsigned int used_size_;

        /// Maximal number of elements on patch.
        unsigned int max_n_elem_;

        /// Distinguishes using of PatchFEValues for storing data of elements or sides.
        MeshObjectType object_type_;
    };

    /// Sub objects of dimensions 1,2,3
    std::array<DimPatchFEValues, 3> dim_fe_vals_;

    uint n_columns_;  ///< Number of columns

};


template <class ValueType>
ValueType ElQ<ValueType>::operator()(FMT_UNUSED const BulkPoint &point) {
    return 0.0;
}

template <>
inline Vector ElQ<Vector>::operator()(FMT_UNUSED const BulkPoint &point) {
    Vector vect; vect.zeros();
    return vect;
}

template <>
inline Tensor ElQ<Tensor>::operator()(FMT_UNUSED const BulkPoint &point) {
	Tensor tens; tens.zeros();
    return tens;
}

template <class ValueType>
ValueType FeQ<ValueType>::operator()(FMT_UNUSED unsigned int shape_idx, FMT_UNUSED const BulkPoint &point) {
    return 0.0;
}

template <>
inline Vector FeQ<Vector>::operator()(FMT_UNUSED unsigned int shape_idx, FMT_UNUSED const BulkPoint &point) {
    Vector vect; vect.zeros();
    return vect;
}

template <>
inline Tensor FeQ<Tensor>::operator()(FMT_UNUSED unsigned int shape_idx, FMT_UNUSED const BulkPoint &point) {
	Tensor tens; tens.zeros();
    return tens;
}


#endif /* PATCH_FE_VALUES_HH_ */
