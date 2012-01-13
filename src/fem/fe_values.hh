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

using namespace std;


/**
 * Class FEValuesData holds the arrays of data computed by
 * Mapping and FiniteElement.
 */
template<unsigned int dim, unsigned int spacedim>
class FEValuesData
{
public:

    /**
     * Resize the data arrays.
     */
    void allocate(unsigned int size, UpdateFlags flags, bool is_scalar = true);



    /**
     * Values at quadrature points of the product of the Jacobian
     * determinant of the mapping and the weight at the particular
     * quadrature point.
     */
    vector<double> JxW_values;

    /**
     * Jacobians of the mapping at the quadrature points.
     */
    vector<arma::mat::fixed<spacedim,dim> > jacobians;

    /**
     * Determinants of Jacobians at quadrature points.
     */
    vector<double> determinants;

    /**
     * Inverse Jacobians at the quadrature points.
     */
    vector<arma::mat::fixed<dim,spacedim> > inverse_jacobians;

    /**
     * Coordinates of quadrature points in the actual cell coordinate system.
     */
    vector<arma::vec::fixed<spacedim> > points;

    /**
     * Shape functions evaluated at the quadrature points.
     */
    vector<arma::vec> shape_values;

    /**
     * Gradients of shape functions evaluated at the quadrature points.
     * Each row of the matrix contains the gradient of one shape function.
     */
    vector<arma::mat> shape_gradients;

    /**
     * Shape functions (for vectorial finite elements) evaluated at
     * quadrature points.
     */
    vector<vector<arma::vec::fixed<spacedim> > > shape_vectors;

    /**
     * Gradients of shape functions (for vectorial finite elements).
     */
    vector<vector<arma::mat::fixed<spacedim,spacedim> > > shape_grad_vectors;

    /**
     * Normal vectors to the element at the quadrature points lying
     * on a side.
     */
    vector<arma::vec::fixed<spacedim> > normal_vectors;

    /**
     * Flags that indicate which finite element quantities are to be computed.
     */
    UpdateFlags update_flags;

    /**
    * Iterator to the last reinit-ed cell.
    */
    typename DOFHandler<dim,spacedim>::CellIterator *present_cell;

};

/**
* Base class for FEValues and FESideValues
*/
template<unsigned int dim, unsigned int spacedim>
class FEValuesBase
{
public:

    /**
     * Default constructor
     */
    FEValuesBase();

    /**
     * Allocates space for computed data.
     */
    void allocate(Mapping<dim,spacedim> &_mapping,
            Quadrature<dim> &_quadrature,
            FiniteElement<dim,spacedim> &_fe,
            UpdateFlags flags);
    
    /**
    * Determine quantities to be recomputed on each cell.
    */
    UpdateFlags update_each(UpdateFlags flags);

    /**
     * Return the value of the @p function_no-th shape function at
     * the @p point_no-th quadrature point.
     */
    const double shape_value(const unsigned int function_no, const unsigned int point_no);

    /**
     * Return the gradient of the @p function_no-th shape function at
     * the @p point_no-th quadrature point.
     */
    const arma::vec::fixed<spacedim> shape_grad(const unsigned int function_no, const unsigned int point_no);

    /**
     * For vectorial finite elements.
     * Return the value of the @p function_no-th shape function at
     * the @p point_no-th quadrature point.
     */
    const arma::vec::fixed<spacedim> shape_vector(const unsigned int function_no, const unsigned int point_no);

    /**
     * For vectorial finite elements.
     * Return the gradient of the @p function_no-th shape function at
     * the @p point_no-th quadrature point.
     */
    const arma::mat::fixed<spacedim,spacedim> shape_grad_vector(const unsigned int function_no, const unsigned int point_no);

    /**
     * Return the relative volume change of the cell (Jacobian determinant).
     * If dim==spacedim then the sign may be negative, otherwise the
     * result is a positive number.
     */
    const double determinant(const unsigned int point_no);

    /**
     * Return the product of Jacobian determinant and the quadrature
     * weight at given quadrature point.
     */
    const double JxW(const unsigned int point_no);

    /**
     * Return coordinates of the quadrature point in the actual cell system.
     */
    const arma::vec::fixed<spacedim> point(const unsigned int point_no);

    /**
     * Returns the normal vector to a side at given quadrature point.
     */
    const arma::vec::fixed<spacedim> normal_vector(unsigned int point_no);


protected:

    /**
     * The mapping from the reference cell to the actual cell.
     */
    Mapping<dim,spacedim> *mapping;

    /**
     * The quadrature rule used to calculate integrals.
     */
    Quadrature<dim> *quadrature;

    /**
     * The used finite element.
     */
    FiniteElement<dim,spacedim> *fe;
    
    /**
     * Precomputed mapping data.
     */
    MappingInternalData *mapping_data;

    /**
     * Precomputed finite element data.
     */
    FEInternalData *fe_data;

    /**
     * Data computed by the mapping and finite element.
     */
    FEValuesData<dim,spacedim> data;
};




/**
 * FEValues takes care of the calculation of finite element data on
 * the actual cell such as values of shape functions at quadrature
 * points, gradients of shape functions, Jacobians of the mapping
 * from the reference cell etc.
 * @param dim      dimension of the reference cell
 * @param spacedim dimension of the Euclidean space where the actual
 *                 cell lives
 *
 * TODO: Implement UpdateFlags to indicate which data should be
 * recomputed at each cell.
 */
template<unsigned int dim, unsigned int spacedim>
class FEValues : public FEValuesBase<dim,spacedim>
{
public:

    FEValues(Mapping<dim,spacedim> &_mapping,
             Quadrature<dim> &_quadrature,
             FiniteElement<dim,spacedim> &_fe,
             UpdateFlags _flags);

    /**
     * Update cell-dependent data (gradients, Jacobians etc.)
     */
    void reinit(typename DOFHandler<dim,spacedim>::CellIterator &cell);


};




/**
 * FESideValues takes care of the calculation of finite element data on
 * a side of the actual cell such as values of shape functions at quadrature
 * points, gradients of shape functions, Jacobians of the mapping
 * from the reference cell etc.
 * @param dim      dimension of the reference cell
 * @param spacedim dimension of the Euclidean space where the actual
 *                 cell lives
 *
 * TODO: Implement UpdateFlags to indicate which data should be
 * recomputed at each cell.
 */
template<unsigned int dim, unsigned int spacedim>
class FESideValues : public FEValuesBase<dim,spacedim>
{

public:

    /**
     * Constructor. Initializes structures and calculates
     * cell-independent data.
     */
    FESideValues(Mapping<dim,spacedim> &_mapping,
             Quadrature<dim-1> &_sub_quadrature,
             FiniteElement<dim,spacedim> &_fe,
             UpdateFlags flags);

    ~FESideValues();

    /**
     * Update cell-dependent data (gradients, Jacobians etc.)
     */
    void reinit(typename DOFHandler<dim,spacedim>::CellIterator &cell,
                Side *side);


private:

    /**
     * Quadrature for the integration on the element sides.
     */
    const Quadrature<dim-1> *sub_quadrature;

};



#endif /* FE_VALUES_HH_ */
