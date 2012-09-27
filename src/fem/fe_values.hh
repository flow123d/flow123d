/*!
 *
 * Copyright (C) 2007 Technical University of Liberec.  All rights reserved.
 *
 * Please make a following refer to Flow123d on your project site if you use the program for any purpose,
 * especially for academic research:
 * Flow123d, Research Centre: Advanced Remedial Technologies, Technical University of Liberec, Czech Republic
 *
 * This program is free software; you can redistribute it and/or modify it under the terms
 * of the GNU General Public License version 3 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 021110-1307, USA.
 *
 *
 * $Id$
 * $Revision$
 * $LastChangedBy$
 * $LastChangedDate$
 *
 * @file
 * @brief Class FEValues calculates finite element data on the actual
 *        cells such as shape function values, gradients, Jacobian of
 *        the mapping from the reference cell etc.
 * @author Jan Stebel
 */

#ifndef FE_VALUES_HH_
#define FE_VALUES_HH_

#include <armadillo>
#include <vector>
#include "fem/update_flags.hh"

template<unsigned int dim, unsigned int spacedim> class DOFHandler;
template<unsigned int dim> class Quadrature;
template<unsigned int dim, unsigned int spacedim> class FEValuesBase;
template<unsigned int dim, unsigned int spacedim> class Mapping;

struct FEInternalData;
struct MappingInternalData;
class Side;



/**
 * @brief Class FEValuesData holds the arrays of data computed by
 * Mapping and FiniteElement.
 */
template<unsigned int dim, unsigned int spacedim>
class FEValuesData
{
public:

    /**
     * @brief Resize the data arrays.
     * @param size Number of quadrature points.
     * @param flags Update flags to be stores.
     * @param is_scalar If true, the structures for scalar values are allocated. Otherwise the vectorial structures are used.
     */
    void allocate(unsigned int size, UpdateFlags flags, bool is_scalar = true);



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
     * @brief Shape functions evaluated at the quadrature points.
     */
    std::vector<arma::vec> shape_values;

    /**
     * @brief Gradients of shape functions evaluated at the quadrature points.
     *
     * Each row of the matrix contains the gradient of one shape function.
     */
    std::vector<arma::mat> shape_gradients;

    /**
     * @brief Shape functions (for vectorial finite elements) evaluated at
     * quadrature points.
     */
    std::vector<std::vector<arma::vec::fixed<spacedim> > > shape_vectors;

    /**
     * @brief Gradients of shape functions (for vectorial finite elements).
     */
    std::vector<std::vector<arma::mat::fixed<spacedim,spacedim> > > shape_grad_vectors;

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
    typename DOFHandler<dim,spacedim>::CellIterator *present_cell;

};


/**
 * @brief Abstract base class with certain methods independent of the template parameter @p dim.
 */
template<unsigned int spacedim>
class FEValuesSpaceBase
{
public:
    /**
     * @brief Return the value of the @p function_no-th shape function at
     * the @p point_no-th quadrature point.
     *
     * @param function_no Number of the shape function.
     * @param point_no Number of the quadrature point.
     */
    virtual const double shape_value(const unsigned int function_no, const unsigned int point_no) = 0;

    /**
     * @brief Return the gradient of the @p function_no-th shape function at
     * the @p point_no-th quadrature point.
     *
     * @param function_no Number of the shape function.
     * @param point_no Number of the quadrature point.
     */
    virtual const arma::vec::fixed<spacedim> shape_grad(const unsigned int function_no, const unsigned int point_no) = 0;

    /**
     * @brief Return the product of Jacobian determinant and the quadrature
     * weight at given quadrature point.
     *
     * @param point_no Number of the quadrature point.
     */
    virtual const double JxW(const unsigned int point_no) = 0;

    /**
     * @brief Returns the normal vector to a side at given quadrature point.
     *
     * @param point_no Number of the quadrature point.
     */
    virtual const arma::vec::fixed<spacedim> normal_vector(unsigned int point_no) = 0;

    /**
     * @brief Returns the number of shape functions.
     */
    virtual const unsigned int n_dofs() = 0;

};

/**
 * @brief Base class for FEValues and FESideValues
 */
template<unsigned int dim, unsigned int spacedim>
class FEValuesBase : public FEValuesSpaceBase<spacedim>
{
public:

    /**
     * @brief Default constructor
     */
    FEValuesBase();


    /**
     * Correct deallocation of objects created by 'initialize' methods.
     */
    virtual ~FEValuesBase();


    /**
     * @brief Allocates space for computed data.
     *
     * @param _mapping The mapping between reference and actual cell.
     * @param _quadrature The quadrature rule.
     * @param _fe The finite element.
     * @param flags The update flags.
     */
    void allocate(Mapping<dim,spacedim> &_mapping,
            Quadrature<dim> &_quadrature,
            FiniteElement<dim,spacedim> &_fe,
            UpdateFlags flags);
    
    /**
     * @brief Determine quantities to be recomputed on each cell.
     *
     * @param flags Flags that indicate what has to be recomputed.
     */
    UpdateFlags update_each(UpdateFlags flags);

    /**
     * @brief Return the value of the @p function_no-th shape function at
     * the @p point_no-th quadrature point.
     *
     * @param function_no Number of the shape function.
     * @param point_no Number of the quadrature point.
     */
    const double shape_value(const unsigned int function_no, const unsigned int point_no);

    /**
     * @brief Return the gradient of the @p function_no-th shape function at
     * the @p point_no-th quadrature point.
     *
     * @param function_no Number of the shape function.
     * @param point_no Number of the quadrature point.
     */
    const arma::vec::fixed<spacedim> shape_grad(const unsigned int function_no, const unsigned int point_no);

    /**
     * @brief Return the value of the @p function_no-th shape function at
     * the @p point_no-th quadrature point.
     *
     * For vectorial finite elements.
     *
     * @param function_no Number of the shape function.
     * @param point_no Number of the quadrature point.
     */
    const arma::vec::fixed<spacedim> shape_vector(const unsigned int function_no, const unsigned int point_no);

    /**
     * @brief Return the gradient of the @p function_no-th shape function at
     * the @p point_no-th quadrature point.
     *
     * For vectorial finite elements.
     *
     * @param function_no Number of the shape function.
     * @param point_no Number of the quadrature point.
     */
    const arma::mat::fixed<spacedim,spacedim> shape_grad_vector(const unsigned int function_no, const unsigned int point_no);

    /**
     * @brief Return the relative volume change of the cell (Jacobian determinant).
     *
     * If dim==spacedim then the sign may be negative, otherwise the
     * result is a positive number.
     *
     * @param point_no Number of the quadrature point.
     */
    const double determinant(const unsigned int point_no);

    /**
     * @brief Return the product of Jacobian determinant and the quadrature
     * weight at given quadrature point.
     *
     * @param point_no Number of the quadrature point.
     */
    const double JxW(const unsigned int point_no);

    /**
     * @brief Return coordinates of the quadrature point in the actual cell system.
     *
     * @param point_no Number of the quadrature point.
     */
    const arma::vec::fixed<spacedim> point(const unsigned int point_no);

    /**
     * @brief Returns the normal vector to a side at given quadrature point.
     *
     * @param point_no Number of the quadrature point.
     */
    const arma::vec::fixed<spacedim> normal_vector(unsigned int point_no);

    /**
     * @brief Returns the number of quadrature points.
     */
    const unsigned int n_points();

    /**
     * @brief Returns the number of shape functions.
     */
    const unsigned int n_dofs();


protected:

    /**
     * @brief The mapping from the reference cell to the actual cell.
     */
    Mapping<dim,spacedim> *mapping;

    /**
     * @brief The quadrature rule used to calculate integrals.
     */
    Quadrature<dim> *quadrature;

    /**
     * @brief The used finite element.
     */
    FiniteElement<dim,spacedim> *fe;
    
    /**
     * @brief Precomputed mapping data.
     */
    MappingInternalData *mapping_data;

    /**
     * @brief Precomputed finite element data.
     */
    FEInternalData *fe_data;

    /**
     * @brief Data computed by the mapping and finite element.
     */
    FEValuesData<dim,spacedim> data;
};




/**
 * @brief Calculates finite element data on the actual cell.
 *
 * FEValues takes care of the calculation of finite element data on
 * the actual cell such as values of shape functions at quadrature
 * points, gradients of shape functions, Jacobians of the mapping
 * from the reference cell etc.
 * @param dim      Dimension of the reference cell.
 * @param spacedim Dimension of the Euclidean space where the actual
 *                 cell lives.
 */
template<unsigned int dim, unsigned int spacedim>
class FEValues : public FEValuesBase<dim,spacedim>
{
public:

	/**
	 * @brief Constructor.
	 *
	 * Initializes structures and calculates
     * cell-independent data.
	 *
	 * @param _mapping The mapping between the reference and actual cell.
	 * @param _quadrature The quadrature rule.
	 * @param _fe The finite element.
	 * @param _flags The update flags.
	 */
    FEValues(Mapping<dim,spacedim> &_mapping,
             Quadrature<dim> &_quadrature,
             FiniteElement<dim,spacedim> &_fe,
             UpdateFlags _flags);

    /**
     * @brief Update cell-dependent data (gradients, Jacobians etc.)
     *
     * @param cell The actual cell.
     */
    void reinit(typename DOFHandler<dim,spacedim>::CellIterator &cell);


};




/**
 * @brief Calculates finite element data on a side.
 *
 * FESideValues takes care of the calculation of finite element data on
 * a side of the actual cell such as values of shape functions at quadrature
 * points, gradients of shape functions, Jacobians of the mapping
 * from the reference cell etc.
 * @param dim      Dimension of the reference cell.
 * @param spacedim Dimension of the Euclidean space where the actual
 *                 cell lives.
 */
template<unsigned int dim, unsigned int spacedim>
class FESideValues : public FEValuesBase<dim,spacedim>
{

public:

    /**
     * @brief Constructor.
     *
     * Initializes structures and calculates
     * cell-independent data.
     *
     * @param _mapping The mapping between the reference and actual cell.
	 * @param _sub_quadrature The quadrature rule on the side.
	 * @param _fe The finite element.
	 * @param _flags The update flags.
     */
    FESideValues(Mapping<dim,spacedim> &_mapping,
             Quadrature<dim-1> &_sub_quadrature,
             FiniteElement<dim,spacedim> &_fe,
             UpdateFlags flags);

    /// Destructor.
    virtual ~FESideValues();

    /**
     * @brief Update cell-dependent data (gradients, Jacobians etc.)
     *
     * @param cell The actual cell.
     * @param side The side of the cell.
     */
    void reinit(typename DOFHandler<dim,spacedim>::CellIterator &cell,
                SideIter side);


private:

    /**
     * @brief Quadrature for the integration on the element sides.
     */
    const Quadrature<dim-1> *sub_quadrature;

};





#endif /* FE_VALUES_HH_ */
