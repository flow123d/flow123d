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
 * @file    fe_value_handler.hh
 * @brief
 */

#ifndef FE_VALUE_HANDLER_HH_
#define FE_VALUE_HANDLER_HH_

#include "system/index_types.hh"
#include "fields/field_values.hh"
#include "fem/finite_element.hh"
#include "mesh/point.hh"
#include "la/vector_mpi.hh"
#include <armadillo>
#include "system/armor.hh"
#include "tools/mixed.hh"

class DOFHandlerMultiDim;
template <int spacedim> class ElementAccessor;


/// Initialization structure of FEValueHandler class.
struct FEValueInitData
{
	/// DOF handler object
    std::shared_ptr<DOFHandlerMultiDim> dh;
    /// Store data of Field
    VectorMPI data_vec;
    /// number of dofs
    unsigned int ndofs;
    /// number of components
    unsigned int n_comp;
    /// Holds begin of component range of evaluation.
    unsigned int range_begin;
    /// Holds end of component range of evaluation.
    unsigned int range_end;
    /// FiniteElement objects of all dimensions.
    MixedPtr<FiniteElement> mixed_fe;
};

/**
 * Helper class that allows compute finite element values specified by element dimension.
 */
template <int elemdim, int spacedim, class Value>
class FEValueHandler
{
public:
	typedef typename Space<spacedim>::Point Point;

	/// Constructor.
	FEValueHandler();

	/// Initialize data members
	void initialize(FEValueInitData init_data);
    /// Returns one value in one given point.
    typename Value::return_type const &value(const Point &p, const ElementAccessor<spacedim> &elm);
    /// Returns std::vector of scalar values in several points at once.
    void value_list (const Armor::array &point_list, const ElementAccessor<spacedim> &elm,
                       std::vector<typename Value::return_type> &value_list);
    /// Compute real coordinates and weights (use QGauss) for given element
    unsigned int compute_quadrature(std::vector<arma::vec::fixed<3>> & q_points, std::vector<double> & q_weights,
    		const ElementAccessor<spacedim> &elm, unsigned int order=3);

    /// Destructor.
	~FEValueHandler();

private:
	/// DOF handler object
    std::shared_ptr<DOFHandlerMultiDim> dh_;
    /// Store data of Field
    VectorMPI data_vec_;
    /// Last value, prevents passing large values (vectors) by value.
    Value value_;
    typename Value::return_type r_value_;
    /// Begin of dof range of actual component
    unsigned int range_begin_;
    /// End of dof range of actual component
    unsigned int range_end_;
    /// Pointer to FiniteElement object used to computing values
    std::shared_ptr<FiniteElement<elemdim>> fe_;
};


/**
 * Specialization for elements with dim==0.
 */
template <int spacedim, class Value>
class FEValueHandler<0, spacedim, Value>
{
public:
	typedef typename Space<spacedim>::Point Point;

	/// Constructor.
	FEValueHandler()
	: value_(r_value_) {}

	/// Initialize data members
	void initialize(FEValueInitData init_data);
    /// Returns one value in one given point.
    typename Value::return_type const &value(const Point &p, const ElementAccessor<spacedim> &elm) {
    	Armor::array point_list(spacedim, 1, 1);
    	point_list.set(0) = Armor::ArmaVec<double,spacedim>( p );
    	std::vector<typename Value::return_type> v_list;
    	v_list.push_back(r_value_);
    	this->value_list(point_list, elm, v_list);
    	this->r_value_ = v_list[0];
    	return this->r_value_;
    }

    /// Returns std::vector of scalar values in several points at once.
    void value_list (const Armor::array &point_list, const ElementAccessor<spacedim> &elm,
                       std::vector<typename Value::return_type> &value_list);

    /// Destructor.
	~FEValueHandler() {}

private:
	/// DOF handler object
    std::shared_ptr<DOFHandlerMultiDim> dh_;
    /// Store data of Field
    VectorMPI data_vec_;
    /// Last value, prevents passing large values (vectors) by value.
    Value value_;
    typename Value::return_type r_value_;
    /// Begin of dof range of actual component
    unsigned int range_begin_;
    /// End of dof range of actual component
    unsigned int range_end_;
};



#endif /* FE_VALUE_HANDLER_HH_ */
