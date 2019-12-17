/*!
 *
﻿ * Copyright (C) 2019 Technical University of Liberec.  All rights reserved.
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
 * @file    element_values.hh
 * @brief   Class ElementValues calculates data related to transformation
 *          of reference cell to actual cell (Jacobian, inverse Jacobian,
 *          determinant, point coordinates etc.).
 * @author  Jan Stebel
 */

#ifndef ELEMENT_VALUES_HH_
#define ELEMENT_VALUES_HH_

// #include <string.h>                           // for memcpy
// #include <algorithm>                          // for swap
// #include <new>                                // for operator new[]
// #include <string>                             // for operator<<
#include <vector>                             // for vector
#include "mesh/ref_element.hh"                // for RefElement
#include "mesh/accessors.hh"                  // for ElementAccessor
#include "fem/update_flags.hh"                // for UpdateFlags

class Quadrature;




/**
 * @brief Structure for storing the precomputed element data.
 */
class RefElementData
{
public:
    
    /// Resize vectors to size @p np.
    RefElementData(unsigned int np);
    
    /// Barycentric coordinates of quadrature points.
    std::vector<arma::vec> bar_coords;
    
    /// Quadrature weights.
    std::vector<double> weights;
    
    /// Number of quadrature points.
    unsigned int n_points;
};




/**
 * @brief Class ElementData holds the arrays of data computed by
 * Mapping.
 */
template<unsigned int dim, unsigned int spacedim>
class ElementData
{
public:
    
    /**
     * @brief Resize the data arrays.
     * @param size   Number of quadrature points.
     * @param flags  Update flags to be stores.
     */
    void allocate(unsigned int size, UpdateFlags flags);
    
    /// Print calculated data.
    void print();

    /**
     * @brief Transformed quadrature weights.
     *
     * Values at quadrature points of the product of the Jacobian
     * determinant of the mapping and the weight at the particular
     * quadrature point.
     */
    std::vector<double> JxW_values;

    /**
     * @brief Jacobians of the mapping at the quadrature points.
     */
    std::vector<arma::mat::fixed<spacedim,dim> > jacobians;

    /**
     * @brief Determinants of Jacobians at quadrature points.
     */
    std::vector<double> determinants;

    /**
     * @brief Inverse Jacobians at the quadrature points.
     */
    std::vector<arma::mat::fixed<dim,spacedim> > inverse_jacobians;

    /**
     * @brief Coordinates of quadrature points in the actual cell coordinate system.
     */
    std::vector<arma::vec::fixed<spacedim> > points;

    /**
     * @brief Normal vectors to the element at the quadrature points lying
     * on a side.
     */
    std::vector<arma::vec::fixed<spacedim> > normal_vectors;

    /**
     * @brief Flags that indicate which finite element quantities are to be computed.
     */
    UpdateFlags update_flags;

    /**
     * @brief Iterator to the last reinit-ed cell.
     */
    ElementAccessor<3> present_cell;

};



/**
 * @brief Base class for ElementValues and ElemSideValues
 */
template<unsigned int dim, unsigned int spacedim>
class ElementValuesBase
{
public:

    /**
     * @brief Default constructor
     */
    ElementValuesBase();


    /**
     * Correct deallocation of objects created by 'initialize' methods.
     */
    virtual ~ElementValuesBase();


    /**
     * @brief Allocates space for computed data.
     *
     * @param n_points Number of quadrature points.
     * @param flags    The update flags.
     */
    void allocate(unsigned int n_points, UpdateFlags flags);
    
    /**
     * @brief Determine quantities to be recomputed on each cell.
     *
     * @param flags Flags that indicate what has to be recomputed.
     */
    UpdateFlags update_each(UpdateFlags flags);
    
    /**
     * @brief Return the relative volume change of the cell (Jacobian determinant).
     *
     * If dim==spacedim then the sign may be negative, otherwise the
     * result is a positive number.
     *
     * @param point_no Number of the quadrature point.
     */
    inline double determinant(const unsigned int point_no) const
    {
        ASSERT_LT_DBG(point_no, n_points_);
        return data.determinants[point_no];
    }
    
    /// Return Jacobian matrix at point @p point_no.
    inline const arma::mat::fixed<spacedim,dim> &jacobian(const unsigned int point_no) const
    {
        ASSERT_LT_DBG(point_no, n_points_);
        return data.jacobians[point_no];
    }
    
    /// Return inverse Jacobian matrix at point @p point_no.
    inline const arma::mat::fixed<dim,spacedim> &inverse_jacobian(const unsigned int point_no) const
    {
        ASSERT_LT_DBG(point_no, n_points_);
        return data.inverse_jacobians[point_no];
    }

    /**
     * @brief Return the product of Jacobian determinant and the quadrature
     * weight at given quadrature point.
     *
     * @param point_no Number of the quadrature point.
     */
    inline double JxW(const unsigned int point_no) const
    {
        ASSERT_LT_DBG(point_no, n_points_);
        return data.JxW_values[point_no];
    }

    /**
     * @brief Return coordinates of the quadrature point in the actual cell system.
     *
     * @param point_no Number of the quadrature point.
     */
    inline const arma::vec::fixed<spacedim> &point(const unsigned int point_no) const
    {
        ASSERT_LT_DBG(point_no, n_points_);
        return data.points[point_no];
    }

    /**
	 * @brief Return coordinates of all quadrature points in the actual cell system.
	 *
	 */
	inline const vector<arma::vec::fixed<spacedim> > &point_list() const
	{
	    return data.points;
	}


    /**
     * @brief Returns the normal vector to a side at given quadrature point.
     *
     * @param point_no Number of the quadrature point.
     */
	inline arma::vec::fixed<spacedim> normal_vector(unsigned int point_no)
	{
        ASSERT_LT_DBG(point_no, n_points_);
	    return data.normal_vectors[point_no];
	}
	
    /**
     * @brief Returns the number of quadrature points.
     */
    inline unsigned int n_points()
    { return n_points_; }
    
    /// Return cell at which the values were reinited.
    const ElementAccessor<spacedim> &cell() const
    { return data.present_cell; }



protected:
    
    /// Precompute data on reference element.
    RefElementData *init_ref_data(const Quadrature &q);
    
    /** @brief Number of integration points. */
    unsigned int n_points_;

    /**
     * @brief Data computed by the mapping.
     */
    ElementData<dim,spacedim> data;
    
};




/**
 * @brief Calculates element data on the actual cell.
 *
 * @param dim      Dimension of the reference cell.
 * @param spacedim Dimension of the Euclidean space where the actual
 *                 cell lives.
 */
template<unsigned int dim, unsigned int spacedim>
class ElementValues : public ElementValuesBase<dim,spacedim>
{
public:

	/**
	 * @brief Constructor.
	 *
	 * Initializes structures and calculates
     * cell-independent data.
	 *
	 * @param _quadrature The quadrature rule.
	 * @param _flags The update flags.
	 */
    ElementValues(Quadrature &_quadrature,
             UpdateFlags _flags);
    
    ~ElementValues();

    /**
     * @brief Update cell-dependent data (gradients, Jacobians etc.)
     *
     * @param cell The actual cell.
     */
    void reinit(const ElementAccessor<3> &cell);
    
    /// Return quadrature.
    const Quadrature &quadrature() const
    { return *quadrature_; }
    
    
    
private:
    
    /// Compute data from reference cell and using MappingP1.
    void fill_data();
    
    /**
     * @brief The quadrature rule used to calculate integrals.
     */
    Quadrature *quadrature_;
    
    /**
     * @brief Precomputed element data.
     */
    RefElementData *ref_data;


};




/**
 * @brief Calculates element data on a side.
 *
 * @param dim      Dimension of the reference cell.
 * @param spacedim Dimension of the Euclidean space where the actual
 *                 cell lives.
 */
template<unsigned int dim, unsigned int spacedim>
class ElemSideValues : public ElementValuesBase<dim,spacedim>
{

public:

    /**
     * @brief Constructor.
     *
     * Initializes structures and calculates
     * cell-independent data.
     *
     * @param _sub_quadrature The quadrature rule on the side (with dimension dim-1).
     * @param flags The update flags.
     */
    ElemSideValues(Quadrature &_sub_quadrature,
             UpdateFlags flags);

    /// Destructor.
    virtual ~ElemSideValues();

    /**
	 * @brief Update cell-dependent data (Jacobians etc.)
	 *
	 * @param cell The actual cell.
	 * @param sid  Number of the side of the cell.
	 */
    void reinit(const ElementAccessor<3> &cell,
        		unsigned int sid);

    /// Return quadrature for given side and its permutation.
    const Quadrature &quadrature(unsigned int sid, unsigned int pid) const
    { return side_quad[sid][pid]; }


private:
    
    /**
     * @brief Calculates the mapping data on a side of a cell.
     */
    void fill_data();

    /**
     * @brief Quadrature for the integration on the element sides.
     */
    const Quadrature *sub_quadrature;

    /// Side quadratures.
    std::vector<std::vector<Quadrature> > side_quad;

    /// Data on reference element (for each side and its permutation).
    RefElementData *side_ref_data[RefElement<dim>::n_sides][RefElement<dim>::n_side_permutations];
    
    /// Current side on which ElemSideValues was recomputed.
    unsigned int side_idx_;
    
};





#endif /* ELEMENT_VALUES_HH_ */
