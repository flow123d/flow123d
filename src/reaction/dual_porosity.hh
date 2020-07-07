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
 * @file    dual_porosity.hh
 * @brief   Class Dual_por_exchange implements the model of dual porosity.
 * 
 * It can be part of the transport model and it computes the concentrations of substances both in 
 * mobile and immobile zone. This model can also work above the sorption model - the sorbed concentration
 * is then computed both from mobile and immobile concentrations. Linear reactions can be define 
 * also in both zones.
 */

#ifndef DUAL_POROSITY_H_
#define DUAL_POROSITY_H_

#include <boost/exception/info.hpp>   // for operator<<, error_info::error_i...
#include <memory>                     // for shared_ptr
#include <vector>
#include "fields/field.hh"            // for Field
#include "fields/field_values.hh"     // for FieldValue<>::Scalar, FieldValue
#include "fields/field_set.hh"
#include "fields/multi_field.hh"
#include "la/vector_mpi.hh"
#include "fields/equation_output.hh"
#include "input/type_base.hh"         // for Array
#include "input/type_generic.hh"      // for Instance
#include "petscvec.h"                 // for Vec, VecScatter, _p_VecScatter
#include "reaction/reaction_term.hh"  // for ReactionTerm

class Mesh;
namespace Input {
	class Record;
	namespace Type {
		class Record;
	}
}


/// Class representing dual porosity model in transport.
class DualPorosity:  public ReactionTerm
{
public:
  typedef ReactionTerm FactoryBaseType;

  /**
   * Static variable for new input data types input
   */
  static const Input::Type::Record & get_input_type();

  /// DualPorosity data
  class EqData : public FieldSet
  {
  public:

    /// Collect all fields
    EqData();

    MultiField<3, FieldValue<3>::Scalar > diffusion_rate_immobile;   ///< Mass transfer coefficients between mobile and immobile pores.
    Field<3, FieldValue<3>::Scalar > porosity_immobile;    ///< Immobile porosity field.
    
    MultiField<3, FieldValue<3>::Scalar> init_conc_immobile; ///< Initial concentrations in the immobile zone. 

    Field<3, FieldValue<3>::Scalar > porosity; ///< Porosity field.
    
    MultiField<3, FieldValue<3>::Scalar>  conc_immobile;    ///< Calculated concentrations in the immobile zone.

    /// Fields indended for output, i.e. all input fields plus those representing solution.
    EquationOutput output_fields;

  };

  /// Constructor.
  DualPorosity(Mesh &init_mesh, Input::Record in_rec);
   
  ///Destructor.
  ~DualPorosity(void);

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
  
  /**
   * Updates the solution according to the dual porosity model.
   */
  void update_solution(void) override;
  
  /// Main output routine.
  void output_data(void) override;
  
  bool evaluate_time_constraint(double &time_constraint) override;
  
protected:
  /**
   * This method disables to use constructor without parameters.
   */
  DualPorosity();

  /// Resolves construction of following reactions.
  void make_reactions();
  
  /// Sets initial condition from input.
  void set_initial_condition();
  /// Initializes field sets.
  void initialize_fields();
  
  void compute_reaction(const DHCellAccessor& dh_cell) override;
  
  /**
   * Pointer to twodimensional array[substance][elements] containing concentrations either in immobile.
   */
  FieldFEScalarVec conc_immobile_fe;

  /**
   * Equation data - all data fields are in this set.
   */
  EqData data_;

  /**
   * Input data set - fields in this set are read from the input file.
   */
  FieldSet input_data_set_;
  
  std::shared_ptr<ReactionTerm> reaction_mobile;       ///< Reaction running in mobile zone
  std::shared_ptr<ReactionTerm> reaction_immobile;     ///< Reaction running in immobile zone
  
  /// Dual porosity computational scheme tolerance. 
  /** According to this tolerance the analytical solution of dual porosity concentrations or
   * simple forward difference approximation of concentrations is chosen for computation.
   */
  double scheme_tolerance_;
  
private:
  /// Registrar of class to factory
  static const int registrar;

};

#endif  //DUAL_POROSITY_H_
