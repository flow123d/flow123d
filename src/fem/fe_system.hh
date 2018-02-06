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
 * @file    fe_system.hh
 * @brief   Class FESystem for compound finite elements.
 * @author  Jan Stebel
 */

#ifndef FE_SYSTEM_HH_
#define FE_SYSTEM_HH_

#include <vector>
#include <memory>

#include "fem/finite_element.hh"
#include "fem/fe_values.hh"



/**
 * @brief Compound finite element on @p dim dimensional simplex.
 *
 * This type of FE is used for vector-valued functions and for systems of equations.
 */
template <unsigned int dim, unsigned int spacedim>
class FESystem : public FiniteElement<dim,spacedim>
{
public:
  
    /**
     * @brief Constructor. FESystem for vector or tensor created from a scalar FE.
     * @param fe Base finite element class.
     * @param t  Type (vector or tensor).
     */
    FESystem(std::shared_ptr<FiniteElement<dim,spacedim> > fe, FEType t);
    
    /**
     * @brief Constructor. FESystem with @p n components created from a scalar FE.
     * @param fe Base finite element class.
     * @param n  Multiplicity (number of components).
     */
    FESystem(const std::shared_ptr<FiniteElement<dim,spacedim> > &fe, unsigned int n);
    
    /**
     * @brief Constructor. FESystem for mixed elements.
     * @param fe Base finite element classes.
     */
    FESystem(std::vector<std::shared_ptr<FiniteElement<dim,spacedim> > > fe);
    
    std::vector<unsigned int> get_scalar_components() const
    { return scalar_components_; }
    
    std::vector<unsigned int> get_vector_components() const
    { return vector_components_; }

    /**
     * @brief The vector variant of basis_value must be implemented but need not be used.
     */
    double basis_value(const unsigned int i, 
                       const arma::vec::fixed<dim> &p, 
                       const unsigned int comp) const override;

    /**
     * @brief The vector variant of basis_grad must be implemented but need not be used.
     */
    arma::vec::fixed<dim> basis_grad(const unsigned int i, 
                                     const arma::vec::fixed<dim> &p, 
                                     const unsigned int comp) const override;
    
    unsigned int n_components() const override { return n_components_; }

    UpdateFlags update_each(UpdateFlags flags) override;

    virtual ~FESystem();

private:

  /**
   * Auxiliary class that stores the relation of dofs in the FESystem to the base FE class.
   * For each dof in FESystem it provides:
   * - base FE class,
   * - index of the basis function within base FE,
   * - component index of the dof in FESystem.
   */
  struct DofComponentData {

    DofComponentData(unsigned int fei, unsigned int bi, unsigned co)
      : fe_index(fei),
        basis_index(bi),
        component_offset(co)
    {};
    
    /// Index of base FE class in the vector fe_.
    unsigned int fe_index;
    /// Index of basis function in the base FE.
    unsigned int basis_index;
    /// Component index in the FESystem.
    unsigned int component_offset;
  };
  
  /// Initialization of the internal structures from the vector of base FE.
  void initialize();
  
  FEInternalData *initialize(const Quadrature<dim> &q) override;

  void fill_fe_values(
      const Quadrature<dim> &q,
      FEInternalData &data,
      FEValuesData<dim,spacedim> &fv_data) override;
  
  void compute_node_matrix() override;
  
  /// Pointers to base FE objects.
  std::vector<std::shared_ptr<FiniteElement<dim,spacedim> > > fe_;
  
  /// Information about dofs.
  std::vector<DofComponentData> fe_dof_indices_;
  
  std::vector<unsigned int> scalar_components_;
  std::vector<unsigned int> vector_components_;
  
  unsigned int n_components_;
  
  /**
   * Auxiliary vector representing permutation of dofs
   * for compatibility with DOFHandler. In DOFHandler,
   * the nodal dofs are ordered in such a way that first
   * come all the dofs on node 1, then all dofs on node 2 etc.
   * On the other hand, in FESystem we take first all dofs from
   * fe_[0], then all dofs from fe_[1] etc.
   * 
   * TODO: Remove this object by modifying the order of fe_dof_indices_
   * or by changing the order of distribution of dofs in DOFHandler.
   */
  std::vector<unsigned int> dof_basis_;
  
};






#endif /* FE_SYSTEM_HH_ */
