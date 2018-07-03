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
 * @brief   Class XFEValues is XFEM extension of FEValues class.
 * @author  Pavel Exner
 */

#ifndef XFE_VALUES_HH_
#define XFE_VALUES_HH_

// #include <string.h>                           // for memcpy
// #include <algorithm>                          // for swap
// #include <new>                                // for operator new[]
// #include <string>                             // for operator<<
#include <vector>                             // for vector
// #include "mesh/ref_element.hh"                // for RefElement
#include "mesh/mesh_types.hh"                 // for ElementFullIter
#include "fem/update_flags.hh"                // for UpdateFlags
// #include "fem/fe_values_views.hh"             // for FEValuesViews
#include "fem/fe_values.hh"                   // for FEValues

template<unsigned int dim> class Quadrature;
template<unsigned int dim> class FiniteElement;
template<unsigned int dim, unsigned int spacedim> class FEValuesBase;
template<unsigned int dim, unsigned int spacedim> class Mapping;

template<int dim, int spacedim> class GlobalEnrichmentFunc;
template<int dim, int spacedim> class XFEMElementData;
template<int dim, int spacedim> class QXFEM;


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
class XFEValues : public FEValuesBase<dim,spacedim>
{
    typedef typename std::shared_ptr<GlobalEnrichmentFunc<dim,spacedim>> EnrichmentPtr;
    typedef typename Space<spacedim>::Point Point;
    
public:

	/**
	 * @brief Constructor.
	 *
	 * Initializes structures and calculates
         * cell-independent data.
	 *
	 * @param _mapping The mapping between the reference and actual cell.
	 * @param _fe The regular finite element.
         * @param pu The finite element representing partition of unity.
	 * @param _flags The update flags.
	 */
    XFEValues(Mapping<dim,spacedim> &_mapping,
              FiniteElement<dim> &_fe,
              FiniteElement<dim> &pu,
              UpdateFlags _flags);

    /**
     * @brief Update cell-dependent data (gradients, Jacobians etc.)
     *
     * @param cell The actual cell.
     */
    void reinit(ElementAccessor<3> &ele,
                XFEMElementData<dim,spacedim> &xdata,
                Quadrature<dim> &_quadrature);
    
    /**
     * @brief Return the relative volume change of the cell (Jacobian determinant).
     *
     * If dim==spacedim then the sign may be negative, otherwise the
     * result is a positive number.
     *
     * @param point_no Number of the quadrature point.
     */
//     inline double determinant(const unsigned int point_no) override
//     {
// TODO:optimize for affine mapping - is constant
//         ASSERT_LT_DBG(point_no, quadrature->size());
//         return data.determinants[point_no];
//     }

    
    /**
     * @brief Returns the normal vector to a side at given quadrature point.
     *
     * @param point_no Number of the quadrature point.
     */
    inline arma::vec::fixed<spacedim> normal_vector(unsigned int point_no) override
    {
        ASSERT_DBG(0).error("Not available in XFEM XFEValues.");
        return arma::vec::fixed<spacedim>();
//         ASSERT_LT_DBG(point_no, quadrature->size());
//         return data.normal_vectors[point_no];
    }

    inline unsigned int n_regular_dofs()
    {   return n_regular_dofs_;}
    
    inline unsigned int n_enriched_dofs()
    {   return n_enriched_dofs_;}

private:
    //enable access to (templated) base class members
    using FEValuesBase<dim,spacedim>::mapping;
    using FEValuesBase<dim,spacedim>::mapping_data;
    using FEValuesBase<dim,spacedim>::fe;
    using FEValuesBase<dim,spacedim>::data;
    using FEValuesBase<dim,spacedim>::quadrature;
    
    void fill_data(const FEInternalData &fe_data) override;
//     void fill_xfem();
    void fill_scalar_xfem_single();
    void fill_vec_piola_xfem_single();
    
    /// Awful HACK function for getting quadrature based on Singularity0D or Singularity1D
    std::shared_ptr<QXFEM<dim,spacedim>> qxfem_side(ElementAccessor<3> &ele, unsigned int sid);
    
    /// Element data cache for SGFEM interpolation.
    /** Keeps integrals of enrichment functions over sides.
     * This optimization is necessary when outputing on refined output meshes.
     */
    struct EleCache{
        typedef std::vector<std::vector<double>> EnrDofValues;
        std::map<unsigned int, EnrDofValues> enr_dof_values;
    };
    static EleCache ele_cache;
    
    /// Finite element used as Partition of Unity.
    FiniteElement<dim> *pu;
    /// XFEM Element data.
    XFEMElementData<dim,spacedim> *xdata;
    
    unsigned int n_regular_dofs_,   ///< Number of regular dofs of the used FE.
                 n_enriched_dofs_;  ///< Number of enriched dofs, given by FE, enrichment technique and number of enrichments.
    
    /// Vector of enrichments.
    std::vector<EnrichmentPtr> enr;

};






#endif /* FE_VALUES_HH_ */
