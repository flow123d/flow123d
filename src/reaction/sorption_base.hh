/** @brief class Sorption is used to enable simulation of sorption described by either linear or Langmuir isotherm in combination with limited solubility under consideration.
 *
 * Class in this file makes it possible to handle the dataset describing solid phase as either precipitated or sorbed species.
 *
 *
 *
 *
 */
#ifndef SORPTION_BASE_H
#define SORPTION_BASE_H

#include <vector>

#include "fields/field_base.hh"
#include "fields/field_set.hh"
#include "reaction/reaction.hh"

class Isotherm;
class Mesh;
class Distribution;

class SorptionBase:  public ReactionTerm
{
public:
  /**
   *   Static variable for new input data types input
   */
  static Input::Type::Record input_type;

  class EqData : public FieldSet // should be written in class Sorption
  {
  public:
    /**
     * Sorption type specifies a kind of equilibrial description of adsorption.
     */
    static Input::Type::Selection sorption_type_selection;

    /// Collect all fields
    EqData(const string &output_field_name);

    Field<3, FieldValue<3>::EnumVector > sorption_type; ///< Discrete need Selection for initialization.
    Field<3, FieldValue<3>::Scalar > rock_density; ///< Rock matrix density.
    Field<3, FieldValue<3>::Vector > isotherm_mult; ///< Multiplication coefficients (k, omega) for all types of isotherms. Langmuir: c_s = omega * (alpha*c_a)/(1- alpha*c_a), Linear: c_s = k*c_a
    Field<3, FieldValue<3>::Vector > isotherm_other; ///< Langmuir sorption coeficients alpha (in fraction c_s = omega * (alpha*c_a)/(1- alpha*c_a)).
    Field<3, FieldValue<3>::Vector> init_conc_solid; ///< Initial sorbed concentrations. 

    Field<3, FieldValue<3>::Scalar > porosity; ///< Porosity field copied from transport
    
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

  void initialize() override;
  void zero_time_step() override;
  void set_initial_condition();

  /**
   * Prepared to compute sorption inside all of considered elements. 
   * It calls compute_reaction(...) for all the elements controled by concrete processor, when the computation is paralelized.
   */
  virtual void update_solution(void);
  
  static Input::Type::Selection make_output_selection(const string &output_field_name, const string &selection_name)
  {
	  return EqData(output_field_name).output_fields.make_output_field_selection(selection_name)
		.close();
  }

  /**
   * Initialization routines that are done in constructors of descendants.
   * Method data() which access EqData is pure virtual and cannot be called from the base constructor.
   */
  //void data_initialization(void);
  
  /**
   * Creates interpolation table for isotherms.
   */
  void make_tables(void);
  
  void output_data(void) override;
  void output_vector_gather(void) override;
  
  /**
   * Meaningless inherited method.
   */
  //void set_concentration_vector(Vec &vec) override;
    
protected:
  /**
   * This method disables to use constructor without parameters.
   */
  SorptionBase();
  
  void initialize_substance_ids(const std::vector<string> &names, Input::Record in_rec);
  
  /// Initializes private members of sorption from the input record.
  void init_from_input(Input::Record in_rec) override;
  
  /// Initializes field sets.
  void initialize_fields();

  /** Initializes possible following reactions from input record.
   * It should be called after setting mesh, time_governor, distribution and concentration_matrix
   * if there are some setting methods for reactions called (they are not at the moment, so it could be part of init_from_input).
   */
  void make_reactions(Input::Record in_rec);
  
  /**
   * For simulation of sorption in just one element either inside of MOBILE or IMMOBILE pores.
   */
  double **compute_reaction(double **concentrations, int loc_el);
  /**
   *
   */
  virtual void isotherm_reinit(std::vector<Isotherm> &isotherms, const ElementAccessor<3> &elm) = 0;
  
  
  void allocate_output_mpi(void);
  


  EqData *data_;

  /**
   * Number of regions.
   */
  int nr_of_regions;
  /**
   * Temporary nr_of_points can be computed using step_length. Should be |nr_of_region x nr_of_substances| matrix later.
   */
  int nr_of_points;
  /**
   * Molar masses of dissolved species (substances)
   */
  std::vector<double> molar_masses;
  /**
   * Density of the solvent. 
   *  TODO: Could be done region dependent, easily.
   */
  double solvent_density;
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

  Input::Type::Selection output_selection;

  /** Reaction model that follows the sorption.
   */
  ReactionTerm* reaction;
                  
  ///@name members used in output routines
  //@{
  Vec *vconc_solid; ///< PETSC sorbed concentration vector (parallel).
  Vec *vconc_solid_out; ///< PETSC sorbed concentration vector output (gathered - sequential)
  double **conc_solid_out; ///< sorbed concentration array output (gathered - sequential)  
  //@}
};

#endif  //SORPTION_BASE_H
