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

#include <vector>

#include "fields/field_algo_base.hh"
#include "fields/field_set.hh"
#include "fields/multi_field.hh"
#include "fields/vec_seq_double.hh"
#include "reaction/reaction_term.hh"

class Isotherm;
class Mesh;

class SorptionBase:  public ReactionTerm
{
public:
    TYPEDEF_ERR_INFO( EI_ArrayName, std::string);
    DECLARE_INPUT_EXCEPTION( ExcSubstanceCountMatch, << "The size of the input array " << EI_ArrayName::qval 
                                                     << " does not match the number of substances.");
    
  /**
   *   Static variable for new input data types input
   */
  static const Input::Type::Record & get_input_type();
  
  static Input::Type::Selection make_output_selection(const string &output_field_name, const string &selection_name)
  {
      return EqData(output_field_name).output_fields
        .make_output_field_selection(selection_name, "desc")
        .close();
  }

  class EqData : public FieldSet
  {
  public:
    /**
     * Sorption type specifies a kind of equilibrial description of adsorption.
     */
    static const Input::Type::Selection & get_sorption_type_selection();

    /// Collect all fields
    EqData(const string &output_field_name);

    Field<3, FieldValue<3>::EnumVector > sorption_type; ///< Discrete need Selection for initialization.
    Field<3, FieldValue<3>::Scalar > rock_density;      ///< Rock matrix density.
    
    /// Multiplication coefficients (k, omega) for all types of isotherms. 
    /** Langmuir: c_s = omega * (alpha*c_a)/(1- alpha*c_a), Linear: c_s = k*c_a */
    Field<3, FieldValue<3>::Vector > isotherm_mult;  
    /// Langmuir sorption coeficients alpha (in fraction c_s = omega * (alpha*c_a)/(1- alpha*c_a)).
    Field<3, FieldValue<3>::Vector > isotherm_other; 
    
    Field<3, FieldValue<3>::Vector> init_conc_solid;    ///< Initial sorbed concentrations. 
    Field<3, FieldValue<3>::Scalar > porosity;          ///< Porosity field copied from transport.
    
    MultiField<3, FieldValue<3>::Scalar>  conc_solid;    ///< Calculated sorbed concentrations, for output only.

    /// Input data set - fields in this set are read from the input file.
    FieldSet input_data_set_;

    /// Fields indended for output, i.e. all input fields plus those representing solution.
    FieldSet output_fields;

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

  ///Reads and sets initial condition for concentration in solid.
  void set_initial_condition();
    
    /// Allocates petsc vectors, prepares them for output and creates vector scatter.
  void allocate_output_mpi(void);
  
  /// Gathers all the parallel vectors to enable them to be output.
  void output_vector_gather(void) override;
  
  /**
   * For simulation of sorption in just one element either inside of MOBILE or IMMOBILE pores.
   */
  double **compute_reaction(double **concentrations, int loc_el);
  
  /// Reinitializes the isotherm.
  /**
   * On data change the isotherm is recomputed, possibly new interpolation table is made.
   * Pure virtual method.
   */
  virtual void isotherm_reinit(std::vector<Isotherm> &isotherms, const ElementAccessor<3> &elm) = 0;
  
    /**
   * Creates interpolation table for isotherms.
   */
  void make_tables(void);
  

  /// Pointer to equation data. The object is constructed in descendants.
  EqData *data_;

  /**
   * Temporary nr_of_points can be computed using step_length. Should be |nr_of_region x nr_of_substances| matrix later.
   */
  unsigned int n_interpolation_steps_;
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
   * Three dimensional array contains intersections between isotherms and mass balance lines. 
   * It describes behaviour of sorbents in mobile pores of various rock matrix enviroments.
   * Up to |nr_of_region x nr_of_substances x n_points| doubles. Because of equidistant step 
   * lenght in cocidered system of coordinates, just function values are stored.
   */
  std::vector<std::vector<Isotherm> > isotherms;
  
  unsigned int n_substances_;   //< number of substances that take part in the sorption mode
  
  /// Mapping from local indexing of substances to global.
  std::vector<unsigned int> substance_global_idx_;
  
  /**
   * Array for storage infos about sorbed species concentrations.
   */
  double** conc_solid;
  
  Input::Array output_array;

  //Input::Type::Selection output_selection;

  /**
   * Reaction model that follows the sorption.
   */
  std::shared_ptr<ReactionTerm> reaction_liquid;
  std::shared_ptr<ReactionTerm> reaction_solid;
                  
  ///@name members used in output routines
  //@{
  VecScatter vconc_out_scatter; ///< Output vector scatter.
  // TODO: replace vconc_solid + conc_solid by VecSeqDouble, use the same principle as in 'conc_solid_out'
  Vec *vconc_solid; ///< PETSC sorbed concentration vector (parallel).
  std::vector<VectorSeqDouble> conc_solid_out; ///< sorbed concentration array output (gathered - sequential)
  //@}
};

#endif  //SORPTION_BASE_H
