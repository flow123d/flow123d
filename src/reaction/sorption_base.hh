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
 * @file    sorption_base.hh
 * @brief   Class SorptionBase is abstract class representing model of sorption in transport.
 * 
 * The sorption is described by several types of isotherms - linear, Freundlich or Langmuir. 
 * Limited solubility can be considered.
 *
 * Interpolation tables are used to speed up evaluation of isotherms.
 *
 */

#ifndef SORPTION_BASE_H
#define SORPTION_BASE_H


#include <memory>                       // for shared_ptr
#include <string>                       // for string
#include <vector>
#include "reaction/reaction_term.hh"    // for ReactionTerm
#include "fields/field.hh"              // for Field
#include "fields/field_values.hh"       // for FieldValue<>::Scalar, FieldVa...
#include "fields/field_set.hh"
#include "fields/multi_field.hh"
#include "la/vector_mpi.hh"
#include "fields/equation_output.hh"
#include "input/input_exception.hh"     // for DECLARE_INPUT_EXCEPTION, Exce...
#include "input/type_base.hh"           // for Array
#include "input/type_generic.hh"        // for Instance
#include "petscvec.h"                   // for Vec, VecScatter, _p_VecScatter
//
#include "system/exceptions.hh"         // for operator<<, ExcStream, EI

class Isotherm;
class Mesh;
namespace Input {
	class Record;
	namespace Type {
		class Record;
		class Selection;
	}
}
template <int spacedim> class ElementAccessor;
template<unsigned int dim> class InitConditionAssemblySorp;
template<unsigned int dim> class ReactionAssemblySorp;
template< template<IntDim...> class DimAssembly> class GenericAssembly;





class SorptionBase:  public ReactionTerm
{
public:
    TYPEDEF_ERR_INFO( EI_ArrayName, std::string);
    TYPEDEF_ERR_INFO( EI_Subst, unsigned int);
    DECLARE_INPUT_EXCEPTION( ExcSubstanceCountMatch, << "The size of the input array " << EI_ArrayName::qval 
                                                     << " does not match the number of substances.");
    DECLARE_EXCEPTION( ExcNotPositiveScaling,
            << "Scaling parameter in sorption is not positive. Check the input for rock density and molar mass of " << EI_Subst::val << ". substance.\n" );
    
  /**
   *   Static variable for new input data types input
   */
  static const Input::Type::Record & get_input_type();

  static Input::Type::Instance make_output_type(const string &equation_name, const string &output_field_name, const string &output_field_desc )
  {
      return EqFields(output_field_name, output_field_desc).output_fields.make_output_type(equation_name, "");
  }

  class EqFields : public ReactionTerm::EqFields
  {
  public:
    /**
     * Sorption type specifies a kind of equilibrial description of adsorption.
     */
    static const Input::Type::Selection & get_sorption_type_selection();

    /// Collect all fields
    EqFields(const string &output_field_name, const string &output_field_desc);

    MultiField<3, FieldValue<3>::Enum > sorption_type; ///< Discrete need Selection for initialization.
    Field<3, FieldValue<3>::Scalar > rock_density;      ///< Rock matrix density.
    
    /// Multiplication coefficients (k, omega) for all types of isotherms. 
    /** Langmuir: c_s = omega * (alpha*c_a)/(1- alpha*c_a), Linear: c_s = k*c_a */
    MultiField<3, FieldValue<3>::Scalar > distribution_coefficient;  
    /// Langmuir sorption coeficients alpha (in fraction c_s = omega * (alpha*c_a)/(1- alpha*c_a)).
    MultiField<3, FieldValue<3>::Scalar > isotherm_other; 
    
    MultiField<3, FieldValue<3>::Scalar> init_conc_solid;    ///< Initial sorbed concentrations. 
    Field<3, FieldValue<3>::Scalar > porosity;          ///< Porosity field copied from transport.
    
    MultiField<3, FieldValue<3>::Scalar>  conc_solid;   ///< Calculated sorbed concentrations, for output only.
    FieldFEScalarVec conc_solid_fe;                     ///< Underlaying FieldFE for each substance of conc_solid.

    /// Input data set - fields in this set are read from the input file.
    FieldSet input_field_set_;

    /// Fields indended for output, i.e. all input fields plus those representing solution.
    EquationOutput output_fields;

    /// Instances of FieldModel used in assembly methods
    Field<3, FieldValue<3>::Scalar > scale_aqua;
    Field<3, FieldValue<3>::Scalar > scale_sorbed;
    Field<3, FieldValue<3>::Scalar > no_sorbing_surface_cond;

  };


  /// DualPorosity data
  class EqData : public ReactionTerm::EqData
  {
  public:

    /// Collect all fields
    EqData();

    /// Returns pair { quad_order_asm, quad_order_fields}
    inline std::vector<unsigned int> quad_order() const {
        return {0, 0};
    }

    /// Mapping from local indexing of substances to global.
    std::vector<unsigned int> substance_global_idx_;

    unsigned int n_substances_;   ///< number of substances that take part in the sorption mode
    /**
     * Density of the solvent.
     *  TODO: Could be done region dependent, easily.
     */
    double solvent_density_;
    /**
     * Critical concentrations of species dissolved in water.
     */
    std::vector<double> solubility_vec_;
    /**
     * Concentration table limits of species dissolved in water.
     */
    std::vector<double> table_limit_;
    /**
     * Maximum concentration per region.
     * It is used for optimization of interpolation table.
     */
    std::vector<std::vector<double>> max_conc;
    /**
     * Three dimensional array contains intersections between isotherms and mass balance lines.
     * It describes behaviour of sorbents in mobile pores of various rock matrix enviroments.
     * Up to |nr_of_region x nr_of_substances x n_points| doubles. Because of equidistant step
     * lenght in cocidered system of coordinates, just function values are stored.
     */
    std::vector<std::vector<Isotherm> > isotherms;
  };


  /**
   *  Constructor with parameter for initialization of a new declared class member
   */
  SorptionBase(Mesh &init_mesh, Input::Record in_rec);
  /**
   * Destructor.
   */
  virtual ~SorptionBase(void);

  /// Prepares the object to usage.
  /**
   * Allocating memory, reading input, initialization of fields.
   */
  void initialize() override;
  
  /**
   * Does first computation after initialization process.
   * The time is set and initial condition is set and output.
   */
  void zero_time_step() override;

  /// Updates the solution. 
  /**
   * Goes through local distribution of elements and calls @p compute_reaction.
   */
  void update_solution(void) override;
  
  void output_data(void) override;
  
    
protected:
  /**
   * This method disables to use constructor without parameters.
   */
  SorptionBase();
  
  /** Initializes possible following reactions from input record.
   * It should be called after setting mesh, time_governor, distribution and concentration_matrix
   * if there are some setting methods for reactions called (they are not at the moment, so it could be part of init_from_input).
   */
  void make_reactions();
  
  /// Reads names of substances from input and creates indexing to global vector of substance.
  /** Also creates the local vector of molar masses. */
  void initialize_substance_ids();
  
  /// Initializes private members of sorption from the input record.
  void initialize_from_input();
  
  /// Initializes field sets.
  void initialize_fields();

  /// Compute reaction on a single element.
  void compute_reaction(const DHCellAccessor& dh_cell) override;
  
//  /// Reinitializes the isotherm.
//  /**
//   * On data change the isotherm is recomputed, possibly new interpolation table is made.
//   */
//  void isotherm_reinit(unsigned int i_subst, const ElementAccessor<3> &elm);
//
//  /// Calls @p isotherm_reinit for all isotherms.
//  void isotherm_reinit_all(const ElementAccessor<3> &elm);
  
    /**
   * Creates interpolation table for isotherms.
   */
  void make_tables(void);
  
  /// Computes maximal aqueous concentration at the current step.
  void update_max_conc();
  
  /// Sets max conc to zeros on all regins.
  void clear_max_conc();

  /// Initialize FieldModels, method is implemented in descendants.
  virtual void init_field_models() = 0;

  std::shared_ptr<EqFields> eq_fields_;  ///< Pointer to equation fields. The object is constructed in descendants.
  std::shared_ptr<EqData> eq_data_;      ///< Equation data

  /**
   * Temporary nr_of_points can be computed using step_length. Should be |nr_of_region x nr_of_substances| matrix later.
   */
  unsigned int n_interpolation_steps_;
  
  /**
   * Reaction model that follows the sorption.
   */
  std::shared_ptr<ReactionTerm> reaction_liquid;
  std::shared_ptr<ReactionTerm> reaction_solid;
  
  /// general assembly objects, hold assembly objects of appropriate dimension
  GenericAssembly< InitConditionAssemblySorp > * init_condition_assembly_;
  GenericAssembly< ReactionAssemblySorp > * reaction_assembly_;


};

#endif  //SORPTION_BASE_H
