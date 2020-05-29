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
#include "tools/mixed.hh"


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



class FESystemFunctionSpace: public FunctionSpace
{
public:

	/**
	 * @brief Constructor.
	 */
    FESystemFunctionSpace(const std::vector<std::shared_ptr<FunctionSpace> > &fs_vector);

    double basis_value(unsigned int basis_index,
                       const arma::vec &point,
                       unsigned int comp_index = 0
                       ) const override;
    
    const arma::vec basis_grad(unsigned int basis_index,
                               const arma::vec &point,
                               unsigned int comp_index = 0
                              ) const override;

    unsigned int dim() const override { return dim_; }
    
    const std::vector<DofComponentData> &dof_indices() { return dof_indices_; }
    
    virtual ~FESystemFunctionSpace() {};

private:

    /// Function spaces that are put together.
    std::vector<std::shared_ptr<FunctionSpace> > fs_;
    
    /// Indices of basis functions relative to the sub-spaces.
    std::vector<DofComponentData> dof_indices_;
    
    /// Number of basis functions.
    unsigned int dim_;
    
};



/**
 * @brief Compound finite element on @p dim dimensional simplex.
 *
 * This type of FE is used for vector-valued functions and for systems of equations.
 */
template <unsigned int dim>
class FESystem : public FiniteElement<dim>
{
public:
  
    /**
     * @brief Constructor for FEVectorContravariant and FEVectorPiola.
     * @param fe Base finite element class (must be scalar).
     * @param t  Type (must be FEVectorContravariant or FEVectorPiola).
     */
    FESystem(std::shared_ptr<FiniteElement<dim> > fe, FEType t);
    
    /**
     * @brief Constructor for FEVector, FETensor and FEMixedSystem.
     * If @p t == FEVector then @p n must be the space dimension into which
     * the reference cell will be mapped and that @p fe is scalar.
     * If @p t == FETensor then @p n must be square of the space dimension into which
     * the reference cell will be mapped and that @p fe is scalar.
     * If @p t == FEMixedSystem, then @p n is the number of components.
     * @param fe Base finite element class (must be scalar if @p t is FEVector or FETensor).
     * @param t  Type of FESystem (must be either FEVector, FETensor or FEMixedSystem).
     * @param n  Multiplicity (number of components).
     */
    FESystem(const std::shared_ptr<FiniteElement<dim> > &fe, FEType t, unsigned int n);
    
    /**
     * @brief Constructor. FESystem for mixed elements.
     * @param fe Base finite element classes.
     */
    FESystem(std::vector<std::shared_ptr<FiniteElement<dim> > > fe);
    
    std::vector<unsigned int> get_scalar_components() const
    { return scalar_components_; }
    
    std::vector<unsigned int> get_vector_components() const
    { return vector_components_; }
    
    std::vector<unsigned int> get_tensor_components() const
    { return tensor_components_; }

    UpdateFlags update_each(UpdateFlags flags) override;
    
    /// Get barycentric coordinates of the points on the reference element associated with the dofs.
    /// Used in BDDC for unknown reason.
    virtual std::vector< arma::vec::fixed<dim+1> > dof_points() const;


    const std::vector<std::shared_ptr<FiniteElement<dim> > > &fe() const
    { return fe_; }
    
    /// Return dof indices belonging to given sub-FE.
    std::vector<unsigned int> fe_dofs(unsigned int fe_index);
    

private:

  /// Initialization of the internal structures from the vector of base FE.
  void initialize();

  void compute_node_matrix() override;
  
  /// Pointers to base FE objects.
  std::vector<std::shared_ptr<FiniteElement<dim> > > fe_;
  
  std::vector<unsigned int> scalar_components_;
  std::vector<unsigned int> vector_components_;
  std::vector<unsigned int> tensor_components_;
  
};

template<class... Args>
MixedPtr<FESystem> mixed_fe_system(MixedPtr<FiniteElement> fe, Args&&... args)
{
    return MixedPtr<FESystem>(
      std::make_shared<FESystem<0>>(fe[0_d], std::forward<Args>(args)...),
      std::make_shared<FESystem<1>>(fe[1_d], std::forward<Args>(args)...),
      std::make_shared<FESystem<2>>(fe[2_d], std::forward<Args>(args)...),
      std::make_shared<FESystem<3>>(fe[3_d], std::forward<Args>(args)...)
      );
}




#endif /* FE_SYSTEM_HH_ */
