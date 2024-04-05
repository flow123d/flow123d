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
#include "fem/eigen_tools.hh"
#include "fem/patch_point_values.hh"
#include "fem/mapping_p1.hh"
#include "mesh/ref_element.hh"                // for RefElement
#include "mesh/accessors.hh"
#include "fem/update_flags.hh"                // for UpdateFlags
#include "quadrature/quadrature_lib.hh"
#include "fields/eval_subset.hh"

template<unsigned int spacedim> class PatchFEValues;



typedef typename std::vector< std::array<uint, 3> > DimPointTable;  ///< Holds triplet (dim; bulk/side; idx of point in subtable)


template <class ValueType>
class ElQ {
public:
    /// Forbidden default constructor
    ElQ() = delete;

    /// Constructor
    ElQ(PatchPointValues<3> &patch_point_vals, unsigned int begin)
    : patch_point_vals_(patch_point_vals), begin_(begin) {}

    ValueType operator()(FMT_UNUSED const BulkPoint &point);

    ValueType operator()(FMT_UNUSED const SidePoint &point);

private:
    // attributes:
    PatchPointValues<3> &patch_point_vals_;
    unsigned int begin_;       /// Index of the first component of the bulk Quantity. Size is given by ValueType
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

    ValueType operator()(FMT_UNUSED unsigned int shape_idx, FMT_UNUSED const SidePoint &point);

    // Implementation for EdgePoint, SidePoint, and JoinPoint shoud have a common implementation
    // resolving to side values

private:
    // attributes:
    PatchFEValues<3> *fe_values_;
    unsigned int begin_;       /// Index of the first component of the Quantity. Size is given by ValueType
};


template <class ValueType>
class JoinShapeAccessor {
public:
    /// Default constructor
    JoinShapeAccessor()
    : join_idx_(-1) {}

    /**
     * Constructor
     *
     * @param fe_values  Pointer to PatchFEValues object.
     * @param begin      Index of the first component of the Quantity.
     * @param lower_dim  Dimension of bulk (lower-dim) element.
     * @param join_idx   Index function.
     */
    JoinShapeAccessor(PatchFEValues<3> *fe_values, unsigned int begin, unsigned int begin_side, unsigned int lower_dim, unsigned int join_idx);

    /// Return global index of DOF
    inline unsigned int join_idx() const {
        return join_idx_;
    }

    /// Return local index of DOF (on low / high-dim) - should be private method
    inline unsigned int local_idx() const {
        if (this->is_high_dim()) return (join_idx_ - n_dofs_low_);
        else return join_idx_;
    }

    inline unsigned int n_dofs_low() const {
        return n_dofs_low_;
    }

    inline unsigned int n_dofs_high() const {
        return n_dofs_high_;
    }

    inline unsigned int n_dofs_both() const {
        return n_dofs_high_ + n_dofs_low_;
    }

    inline bool is_high_dim() const {
        return (join_idx_ >= n_dofs_low_);
    }

    /// Iterates to next item.
    inline void inc() {
        join_idx_++;
    }

    /// Comparison of accessors.
    bool operator==(const JoinShapeAccessor<ValueType>& other) const {
    	return (join_idx_ == other.join_idx_);
    }


    ValueType operator()(FMT_UNUSED const BulkPoint &point);

    ValueType operator()(FMT_UNUSED const SidePoint &point);

private:
    // attributes:
    PatchFEValues<3> *fe_values_;
    unsigned int begin_;          ///< Index of the first component of the bulk Quantity. Size is given by ValueType
    unsigned int begin_side_;     ///< Index of the first component of the side Quantity. Size is given by ValueType
    unsigned int n_dofs_high_;
    unsigned int n_dofs_low_;
    unsigned int join_idx_;
};


template<unsigned int Dim>
class BulkValues
{
public:
	/// Constructor
	BulkValues(PatchPointValues<3> &patch_point_vals)
	: patch_point_vals_(patch_point_vals) {
	    ASSERT_EQ(patch_point_vals.dim(), Dim);
	}

	/// Return dimension
	inline unsigned int dim() const {
	    return Dim;
	}

    /**
     * @brief Register the product of Jacobian determinant and the quadrature
     * weight at bulk quadrature points.
     *
     * @param quad Quadrature.
     */
    inline ElQ<Scalar> JxW()
    {
        uint begin = patch_point_vals_.operations_[FeBulk::BulkOps::opJxW].result_row();
        return ElQ<Scalar>(patch_point_vals_, begin);
    }

	/// Create bulk accessor of coords entity
    inline ElQ<Vector> coords()
    {
        uint begin = patch_point_vals_.operations_[FeBulk::BulkOps::opCoords].result_row();
        return ElQ<Vector>(patch_point_vals_, begin);
    }

//    inline ElQ<Tensor> jacobian(std::initializer_list<Quadrature *> quad_list)
//    {}

    /// Create bulk accessor of jac determinant entity
    inline ElQ<Scalar> determinant()
    {
        uint begin = patch_point_vals_.operations_[FeBulk::BulkOps::opJacDet].result_row();
        return ElQ<Scalar>(patch_point_vals_, begin);
    }

private:
    PatchPointValues<3> &patch_point_vals_;
};


template<unsigned int Dim>
class SideValues
{
public:
	/// Constructor
	SideValues(PatchPointValues<3> &patch_point_vals)
	: patch_point_vals_(patch_point_vals) {
	    ASSERT_EQ(patch_point_vals.dim(), Dim);
	}

	/// Return dimension
	inline unsigned int dim() const {
	    return Dim;
	}

    /// Same as BulkValues::JxW but register at side quadrature points.
    inline ElQ<Scalar> JxW()
    {
        uint begin = patch_point_vals_.operations_[FeSide::SideOps::opJxW].result_row();
        return ElQ<Scalar>(patch_point_vals_, begin);
    }

    /**
     * @brief Register the normal vector to a side at side quadrature points.
     *
     * @param quad Quadrature.
     */
	inline ElQ<Vector> normal_vector()
	{
        uint begin = patch_point_vals_.operations_[FeSide::SideOps::opNormalVec].result_row();
        return ElQ<Vector>(patch_point_vals_, begin);
	}

	/// Create side accessor of coords entity
    inline ElQ<Vector> coords()
    {
        uint begin = patch_point_vals_.operations_[FeSide::SideOps::opCoords].result_row();
        return ElQ<Vector>(patch_point_vals_, begin);
    }

    /// Create bulk accessor of jac determinant entity
    inline ElQ<Scalar> determinant()
    {
        uint begin = patch_point_vals_.operations_[FeSide::SideOps::opSdJacDet].result_row();
        return ElQ<Scalar>(patch_point_vals_, begin);
    }

private:
    PatchPointValues<3> &patch_point_vals_;
};


template<unsigned int spacedim = 3>
class PatchFEValues {
private:
    enum MeshObjectType {
        ElementFE = 0,
		SideFE = 1
    };

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
         * !! Temporary function. Use in fill_data.
         * Set shape value @p val of the @p i_point and @p i_func_comp.
         */
        inline void set_shape_value(unsigned int i_point, unsigned int i_func_comp, double val)
        {
            element_data_[patch_data_idx_].shape_values_[i_point][i_func_comp] = val;
        }

        /**
         * !! Temporary function. Use in fill_data.
         * Set shape gradient @p val of the @p i_point and @p i_func_comp.
         */
        inline void set_shape_gradient(unsigned int i_point, unsigned int i_func_comp, arma::vec::fixed<spacedim> val)
        {
            element_data_[patch_data_idx_].shape_gradients_[i_point][i_func_comp] = val;
        }

        /**
         * !! Temporary function. Use in fill_data.
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
            return element_data_[patch_data_idx_].shape_values_[point_no][function_no];
        }

        /**
         * !! Temporary function. Use in fill_data.
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
            return element_data_[patch_data_idx_].shape_gradients_[point_no][function_no];
        }

        /**
         * @brief Return the product of Jacobian determinant and the quadrature
         * weight at given quadrature point.
         *
         * @param p BulkPoint corresponds to the quadrature point.
         */
        inline double JxW(const BulkPoint &p)
        {
            unsigned int patch_data_idx = element_patch_map_.find(p.elem_patch_idx())->second;
            return element_data_[patch_data_idx].elm_values_->JxW(p.eval_point_idx());
        }

        /**
         * @brief Return the product of Jacobian determinant and the quadrature
         * weight at given quadrature point.
         *
         * @param p SidePoint corresponds to the quadrature point.
         */
        inline double JxW(const SidePoint &p)
        {
            unsigned int patch_data_idx = element_patch_map_.find(p.elem_patch_idx())->second + p.side_idx();
            return element_data_[patch_data_idx].elm_values_->side_JxW(p.local_point_idx());
        }

        /**
         * @brief Return the value of the @p function_no-th shape function at
         * the @p p quadrature point.
         *
         * @param function_no Number of the shape function.
         * @param p BulkPoint corresponds to the quadrature point.
         */
        inline double shape_value(const unsigned int function_no, const BulkPoint &p) const
        {
            ASSERT_LT(function_no, this->n_dofs_);
            unsigned int patch_data_idx = element_patch_map_.find(p.elem_patch_idx())->second;
            return element_data_[patch_data_idx].shape_values_[p.eval_point_idx()][function_no];
        }

        /**
         * @brief Return the value of the @p function_no-th shape function at
         * the @p p quadrature point.
         *
         * @param function_no Number of the shape function.
         * @param p SidePoint corresponds to the quadrature point.
         */
        inline double shape_value(const unsigned int function_no, const SidePoint &p) const
        {
            ASSERT_LT(function_no, this->n_dofs_);
            unsigned int patch_data_idx = element_patch_map_.find(p.elem_patch_idx())->second + p.side_idx();
            return element_data_[patch_data_idx].shape_values_[p.local_point_idx()][function_no];
        }

        /**
         * @brief Return the gradient of the @p function_no-th shape function at
         * the @p p quadrature point.
         *
         * @param function_no Number of the shape function.
         * @param p BulkPoint corresponds to the quadrature point.
         */
        inline arma::vec::fixed<spacedim> shape_grad(const unsigned int function_no, const BulkPoint &p) const
    	{
            ASSERT_LT(function_no, this->n_dofs_);
            unsigned int patch_data_idx = element_patch_map_.find(p.elem_patch_idx())->second;
            return element_data_[patch_data_idx].shape_gradients_[p.eval_point_idx()][function_no];;
        }

        /**
         * @brief Return the gradient of the @p function_no-th shape function at
         * the @p p quadrature point.
         *
         * @param function_no Number of the shape function.
         * @param p SidePoint corresponds to the quadrature point.
         */
        inline arma::vec::fixed<spacedim> shape_grad(const unsigned int function_no, const SidePoint &p) const
    	{
            ASSERT_LT(function_no, this->n_dofs_);
            unsigned int patch_data_idx = element_patch_map_.find(p.elem_patch_idx())->second + p.side_idx();
            return element_data_[patch_data_idx].shape_gradients_[p.local_point_idx()][function_no];;
        }

        /**
         * @brief Returns the normal vector to a side at given quadrature point.
         *
         * @param p SidePoint corresponds to the quadrature point.
         */
    	inline arma::vec::fixed<spacedim> normal_vector(const SidePoint &p)
    	{
            unsigned int patch_data_idx = element_patch_map_.find(p.elem_patch_idx())->second + p.side_idx();
            return element_data_[patch_data_idx].elm_values_->normal_vector(p.local_point_idx());
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

    /// Temporary helper class used in step between usage old a new implementation
	class FuncDef {
    public:
    	FuncDef() {}
        FuncDef(DimPatchFEValues *data, string func_name)
        : point_data_(data), func_name_(func_name) {}
        DimPatchFEValues *point_data_;
	    string func_name_;
    };
public:

    PatchFEValues()
    : dim_fe_vals_({DimPatchFEValues(0), DimPatchFEValues(0), DimPatchFEValues(0)}),
      dim_fe_side_vals_({DimPatchFEValues(0), DimPatchFEValues(0), DimPatchFEValues(0)}),
	  patch_point_vals_bulk_{ {FeBulk::PatchPointValues(1, 0, 1), FeBulk::PatchPointValues(2, 0, 1), FeBulk::PatchPointValues(3, 0, 1)} },
	  patch_point_vals_side_{ {FeSide::PatchPointValues(1, 0, 1), FeSide::PatchPointValues(2, 0, 1), FeSide::PatchPointValues(3, 0, 1)} } {
        used_quads_[0] = false; used_quads_[1] = false;
    }

    PatchFEValues(unsigned int n_quad_points, unsigned int quad_order, MixedPtr<FiniteElement> fe)
    : dim_fe_vals_({DimPatchFEValues(n_quad_points), DimPatchFEValues(n_quad_points), DimPatchFEValues(n_quad_points)}),
      dim_fe_side_vals_({DimPatchFEValues(n_quad_points), DimPatchFEValues(n_quad_points), DimPatchFEValues(n_quad_points)}),
      patch_point_vals_bulk_{ {FeBulk::PatchPointValues(1, quad_order, fe[Dim<1>{}]->n_dofs()),
    	                       FeBulk::PatchPointValues(2, quad_order, fe[Dim<2>{}]->n_dofs()),
                               FeBulk::PatchPointValues(3, quad_order, fe[Dim<3>{}]->n_dofs())} },
      patch_point_vals_side_{ {FeSide::PatchPointValues(1, quad_order, fe[Dim<0>{}]->n_dofs()),
                               FeSide::PatchPointValues(2, quad_order, fe[Dim<1>{}]->n_dofs()),
                               FeSide::PatchPointValues(3, quad_order, fe[Dim<2>{}]->n_dofs())} },
      fe_(fe) {
        used_quads_[0] = false; used_quads_[1] = false;
    }


    /// Destructor
    ~PatchFEValues()
    {}

    /// Return bulk or side quadrature of given dimension
    Quadrature *get_quadrature(uint dim, bool is_bulk) const {
        if (is_bulk) return patch_point_vals_bulk_[dim-1].get_quadrature();
        else return patch_point_vals_side_[dim-1].get_quadrature();
    }

    /**
	 * @brief Initialize structures and calculates cell-independent data.
	 *
	 * @param _quadrature The quadrature rule for the cell associated
     *                    to given finite element or for the cell side.
	 * @param _flags The update flags.
	 */
    template<unsigned int DIM>
    void initialize(Quadrature &_quadrature,
                    UpdateFlags _flags)
    {
        if ( _quadrature.dim() == DIM ) {
            dim_fe_vals_[DIM-1].initialize(_quadrature, *fe_[Dim<DIM>{}], _flags);
            used_quads_[0] = true;
            // new data storing
            patch_point_vals_bulk_[DIM-1].initialize(3); // bulk
        } else {
            dim_fe_side_vals_[DIM-1].initialize(_quadrature, *fe_[Dim<DIM>{}], _flags);
            used_quads_[1] = true;
            // new data storing
            patch_point_vals_side_[DIM-1].initialize(4); // side
        }
    }

    /// Reset PatchpointValues structures
    void reset()
    {
        for (unsigned int i=0; i<3; ++i) {
            if (used_quads_[0]) patch_point_vals_bulk_[i].reset();
            if (used_quads_[1]) patch_point_vals_side_[i].reset();
        }
    }

    /// Reinit data - old data storing, temporary
    void reinit(std::array<PatchElementsList, 4> patch_elements)
    {
        for (unsigned int i=0; i<3; ++i) {
            if (used_quads_[0]) dim_fe_vals_[i].reinit(patch_elements[i+1]);
            if (used_quads_[1]) dim_fe_side_vals_[i].reinit(patch_elements[i+1]);
        }
    }

    /// Reinit data.
    void reinit_patch()
    {
        for (unsigned int i=0; i<3; ++i) {
            if (used_quads_[0]) patch_point_vals_bulk_[i].reinit_patch();
            if (used_quads_[1]) patch_point_vals_side_[i].reinit_patch();
        }
    }

    /**
     * @brief Returns the number of shape functions.
     */
    inline unsigned int n_dofs(unsigned int dim, FMT_UNUSED uint component_idx = 0) const
    {
        ASSERT( (dim>0) && (dim<=3) )(dim).error("Invalid dimension!");
        return dim_fe_vals_[dim-1].n_dofs_;
    }

    /// Return BulkValue object of dimension given by template parameter
    template<unsigned int dim>
    BulkValues<dim> bulk_values() {
        return BulkValues<dim>(patch_point_vals_bulk_[dim-1]);
    }

    /// Return SideValue object of dimension given by template parameter
    template<unsigned int dim>
    SideValues<dim> side_values() {
        return SideValues<dim>(patch_point_vals_side_[dim-1]);
    }

    /**
     * @brief Return the value of the @p function_no-th shape function at
     * the @p p quadrature point.
     *
     * @param quad_list List of quadratures.
     * @param function_no Number of the shape function.
     */
    inline FeQ<Scalar> scalar_shape(Quadrature *quad, uint component_idx = 0)
    {
        uint dim = quad->dim();
        uint begin = patch_point_vals_bulk_[dim-1].add_rows(this->n_dofs(dim, component_idx)); // scalar needs one column
        func_map_[begin] = FuncDef( &dim_fe_vals_[dim-1], "shape_value"); // storing to temporary map

        return FeQ<Scalar>(this, begin);
    }

    /// Same as previous but register at side quadrature points.
    inline FeQ<Scalar> scalar_shape_side(Quadrature *quad, uint component_idx = 0)
    {
        uint dim = quad->dim();
       	uint begin = patch_point_vals_side_[dim].add_rows(this->n_dofs(dim+1, component_idx));  // scalar needs one column
        func_map_side_[begin] = FuncDef( &dim_fe_side_vals_[dim], "shape_value");

        return FeQ<Scalar>(this, begin);
    }

//    inline FeQ<Vector> vector_shape(std::initializer_list<Quadrature *> quad_list)
//    {}

//    inline FeQ<Tensor> tensor_shape(std::initializer_list<Quadrature *> quad_list)
//    {}

    inline FeQ<Vector> grad_scalar_shape(Quadrature *quad, uint component_idx=0)
    {
        uint dim = quad->dim();
       	uint begin = patch_point_vals_bulk_[dim-1].add_rows(3 * this->n_dofs(dim, component_idx)); // Vector needs 3 columns column * n_dofs
        func_map_[begin] = FuncDef( &dim_fe_vals_[dim-1], "shape_grad"); // storing to temporary map

        return FeQ<Vector>(this, begin);
    }

    /// Same as previous but register at side quadrature points.
    inline FeQ<Vector> grad_scalar_shape_side(Quadrature *quad, uint component_idx=0)
    {
        uint dim = quad->dim();
       	uint begin = patch_point_vals_side_[dim].add_rows(3 * this->n_dofs(dim+1, component_idx));  // Vector needs 3 columns column * n_dofs
        func_map_side_[begin] = FuncDef( &dim_fe_side_vals_[dim], "shape_grad");

        return FeQ<Vector>(this, begin);
    }

//    inline FeQ<Tensor> grad_vector_shape(std::initializer_list<Quadrature *> quad_list, unsigned int i_comp=0)
//    {}

    inline Range< JoinShapeAccessor<Scalar> > scalar_join_shape(std::initializer_list<Quadrature *> quad_list)
    {
        std::vector<Quadrature *> quad_vec(quad_list);
        uint dim = quad_vec[0]->dim();
        uint begin=-1, begin_side=-1;

        if (quad_vec[0] != nullptr) {
            begin = patch_point_vals_bulk_[dim-1].add_rows(1); // scalar needs one column
            func_map_[begin] = FuncDef( &dim_fe_vals_[dim-1], "scalar_join_shape"); // storing to temporary map
        }
        if (quad_vec[1] != nullptr) {
            begin_side = patch_point_vals_side_[dim-1].add_rows(1); // scalar needs one column
            func_map_side_[begin_side] = FuncDef( &dim_fe_side_vals_[dim], "scalar_join_shape"); // storing to temporary map
        }

        auto bgn_it = make_iter<JoinShapeAccessor<Scalar>>( JoinShapeAccessor<Scalar>(this, begin, begin_side, quad_vec[0]->dim(), 0) );
        unsigned int end_idx = bgn_it->n_dofs_high() + bgn_it->n_dofs_low();
        auto end_it = make_iter<JoinShapeAccessor<Scalar>>( JoinShapeAccessor<Scalar>(this, begin, begin_side, quad_vec[0]->dim(), end_idx) );
        return Range<JoinShapeAccessor<Scalar>>(bgn_it, end_it);
    }

    /// Resize \p dim_point_table_ if actual size is less than new_size and return reference
    inline DimPointTable &dim_point_table(unsigned int new_size) {
        if (dim_point_table_.size() < new_size) dim_point_table_.resize(new_size);
        for (uint i=0; i<new_size; ++i) dim_point_table_[i][0] = 10; // set invalid dim
        return dim_point_table_;
    }

    /** Following methods are used during update of patch. **/

    /// Resize tables of patch_point_vals_
    void resize_tables(std::vector<std::vector<uint> > dim_sizes) {
        ASSERT_EQ(dim_sizes.size(), 2);
        ASSERT_EQ(dim_sizes[0].size(), 3);

        for (uint i=0; i<3; ++i) {
        	patch_point_vals_bulk_[i].resize_tables(dim_sizes[0][i]);
        	patch_point_vals_side_[i].resize_tables(dim_sizes[1][i]);
        }
    }

    /// Register element to patch_point_vals_ table by dimension of element
    uint register_element(DHCellAccessor cell, uint element_patch_idx) {
        arma::mat coords;
        switch (cell.dim()) {
        case 1:
            coords = MappingP1<1,spacedim>::element_map(cell.elm());
            return patch_point_vals_bulk_[0].register_element(coords, element_patch_idx);
            break;
        case 2:
        	coords = MappingP1<2,spacedim>::element_map(cell.elm());
            return patch_point_vals_bulk_[1].register_element(coords, element_patch_idx);
            break;
        case 3:
        	coords = MappingP1<3,spacedim>::element_map(cell.elm());
            return patch_point_vals_bulk_[2].register_element(coords, element_patch_idx);
            break;
        default:
        	ASSERT(false);
        	return 0;
            break;
        }
    }

    /// Register side to patch_point_vals_ table by dimension of side
    uint register_side(DHCellSide cell_side) {
        arma::mat side_coords(spacedim, cell_side.dim());
        for (unsigned int n=0; n<cell_side.dim(); n++)
            for (unsigned int c=0; c<spacedim; c++)
                side_coords(c,n) = (*cell_side.side().node(n))[c];

        arma::mat elm_coords;
        DHCellAccessor cell = cell_side.cell();
        switch (cell.dim()) {
        case 1:
            elm_coords = MappingP1<1,spacedim>::element_map(cell.elm());
            return patch_point_vals_side_[0].register_side(elm_coords, side_coords);
            break;
        case 2:
            elm_coords = MappingP1<2,spacedim>::element_map(cell.elm());
            return patch_point_vals_side_[1].register_side(elm_coords, side_coords);
            break;
        case 3:
            elm_coords = MappingP1<3,spacedim>::element_map(cell.elm());
            return patch_point_vals_side_[2].register_side(elm_coords, side_coords);
            break;
        default:
            ASSERT(false);
            return 0;
            break;
        }
    }

    /// Register bulk point to patch_point_vals_ table by dimension of element
    uint register_bulk_point(DHCellAccessor cell, uint elem_table_row, uint value_patch_idx) {
        return patch_point_vals_bulk_[cell.dim()-1].register_bulk_point(elem_table_row, value_patch_idx, cell.elm_idx());
    }

    /// Register side point to patch_point_vals_ table by dimension of side
    uint register_side_point(DHCellSide cell_side, uint elem_table_row, uint value_patch_idx) {
        return patch_point_vals_side_[cell_side.dim()-1].register_side_point(elem_table_row, value_patch_idx, cell_side.elem_idx(), cell_side.side_idx());
    }

    /// Temporary development method
    void print(bool points, bool ints, bool only_bulk=true) const {
        for (uint i=0; i<3; ++i)
            patch_point_vals_bulk_[i].print(points, ints);
        if (!only_bulk)
            for (uint i=0; i<3; ++i)
                patch_point_vals_side_[i].print(points, ints);
    }

private:
    /// Sub objects of dimensions 1,2,3
    std::array<DimPatchFEValues, 3> dim_fe_vals_;
    std::array<DimPatchFEValues, 3> dim_fe_side_vals_;
    std::array<FeBulk::PatchPointValues<spacedim>, 3> patch_point_vals_bulk_;
    std::array<FeSide::PatchPointValues<spacedim>, 3> patch_point_vals_side_;
    DimPointTable dim_point_table_;

    MixedPtr<FiniteElement> fe_;         ///< Mixed of shared pointers of FiniteElement object

    ///< Temporary helper objects used in step between usage old a new implementation
    bool used_quads_[2];
    std::map<unsigned int, FuncDef> func_map_;
    std::map<unsigned int, FuncDef> func_map_side_;

    template <class ValueType>
    friend class ElQ;
    template <class ValueType>
    friend class FeQ;
    template <class ValueType>
    friend class JoinShapeAccessor;
};


template <class ValueType>
ValueType ElQ<ValueType>::operator()(const BulkPoint &point) {
    unsigned int value_cache_idx = point.elm_cache_map()->element_eval_point(point.elem_patch_idx(), point.eval_point_idx());
    return patch_point_vals_.scalar_val(begin_, value_cache_idx);
}

template <>
inline Vector ElQ<Vector>::operator()(const BulkPoint &point) {
    unsigned int value_cache_idx = point.elm_cache_map()->element_eval_point(point.elem_patch_idx(), point.eval_point_idx());
    return patch_point_vals_.vector_val(begin_, value_cache_idx);
}

template <>
inline Tensor ElQ<Tensor>::operator()(const BulkPoint &point) {
    unsigned int value_cache_idx = point.elm_cache_map()->element_eval_point(point.elem_patch_idx(), point.eval_point_idx());
    return patch_point_vals_.tensor_val(begin_, value_cache_idx);
}

template <class ValueType>
ValueType ElQ<ValueType>::operator()(const SidePoint &point) {
    unsigned int value_cache_idx = point.elm_cache_map()->element_eval_point(point.elem_patch_idx(), point.eval_point_idx());
    return patch_point_vals_.scalar_val(begin_, value_cache_idx);
}

template <>
inline Vector ElQ<Vector>::operator()(const SidePoint &point) {
    unsigned int value_cache_idx = point.elm_cache_map()->element_eval_point(point.elem_patch_idx(), point.eval_point_idx());
    return patch_point_vals_.vector_val(begin_, value_cache_idx);
}

template <>
inline Tensor ElQ<Tensor>::operator()(const SidePoint &point) {
    unsigned int value_cache_idx = point.elm_cache_map()->element_eval_point(point.elem_patch_idx(), point.eval_point_idx());
    return patch_point_vals_.tensor_val(begin_, value_cache_idx);
}

template <class ValueType>
ValueType FeQ<ValueType>::operator()(unsigned int shape_idx, const BulkPoint &point) {
	auto it = fe_values_->func_map_.find(begin_);
    if (it->second.func_name_ == "shape_value") {
        return it->second.point_data_->shape_value(shape_idx, point);
    } else {
        //ASSERT_PERMANENT(false).error("Should not happen.");
        return 0.0;
    }
}

template <>
inline Vector FeQ<Vector>::operator()(unsigned int shape_idx, const BulkPoint &point) {
	auto it = fe_values_->func_map_.find(begin_);
    if (it->second.func_name_ == "shape_grad") {
        return it->second.point_data_->shape_grad(shape_idx, point);
    } else {
        //ASSERT_PERMANENT(false).error("Should not happen.");
        Vector vect; vect.zeros();
        return vect;
    }
}

template <>
inline Tensor FeQ<Tensor>::operator()(FMT_UNUSED unsigned int shape_idx, FMT_UNUSED const BulkPoint &point) {
	Tensor tens; tens.zeros();
    return tens;
}

template <class ValueType>
ValueType FeQ<ValueType>::operator()(unsigned int shape_idx, const SidePoint &point) {
	auto it = fe_values_->func_map_side_.find(begin_);
    if (it->second.func_name_ == "shape_value") {
        return it->second.point_data_->shape_value(shape_idx, point);
    } else {
        //ASSERT_PERMANENT(false).error("Should not happen.");
        return 0.0;
    }
}

template <>
inline Vector FeQ<Vector>::operator()(unsigned int shape_idx, const SidePoint &point) {
	auto it = fe_values_->func_map_side_.find(begin_);
    if (it->second.func_name_ == "shape_grad") {
        return it->second.point_data_->shape_grad(shape_idx, point);
    } else {
        //ASSERT_PERMANENT(false).error("Should not happen.");
        Vector vect; vect.zeros();
        return vect;
    }
}

template <>
inline Tensor FeQ<Tensor>::operator()(FMT_UNUSED unsigned int shape_idx, FMT_UNUSED const SidePoint &point) {
	Tensor tens; tens.zeros();
    return tens;
}


template <class ValueType>
ValueType JoinShapeAccessor<ValueType>::operator()(const BulkPoint &point) {
	auto it = fe_values_->func_map_.find(begin_);
    if (it->second.func_name_ == "scalar_join_shape") {
        if (this->is_high_dim()) return 0.0;
        else return it->second.point_data_->shape_value(this->local_idx(), point);
    } else {
        //ASSERT_PERMANENT(false).error("Should not happen.");
        return 0.0;
    }
}

template <>
inline Vector JoinShapeAccessor<Vector>::operator()(FMT_UNUSED const BulkPoint &point) {
    Vector vect; vect.zeros();
    return vect;
}

template <>
inline Tensor JoinShapeAccessor<Tensor>::operator()(FMT_UNUSED const BulkPoint &point) {
	Tensor tens; tens.zeros();
    return tens;
}

template <class ValueType>
ValueType JoinShapeAccessor<ValueType>::operator()(const SidePoint &point) {
	auto it = fe_values_->func_map_side_.find(begin_side_);
    if (it->second.func_name_ == "scalar_join_shape") {
        if (this->is_high_dim()) return it->second.point_data_->shape_value(this->local_idx(), point);
        else return 0.0;
    } else {
        //ASSERT_PERMANENT(false).error("Should not happen.");
        return 0.0;
    }
}

template <>
inline Vector JoinShapeAccessor<Vector>::operator()(FMT_UNUSED const SidePoint &point) {
    Vector vect; vect.zeros();
    return vect;
}

template <>
inline Tensor JoinShapeAccessor<Tensor>::operator()(FMT_UNUSED const SidePoint &point) {
	Tensor tens; tens.zeros();
    return tens;
}


#endif /* PATCH_FE_VALUES_HH_ */
