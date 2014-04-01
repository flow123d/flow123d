/** @brief class Sorption is used to enable simulation of sorption described by either linear or Langmuir isotherm in combination with limited solubility under consideration.
 *
 * Class in this file makes it possible to handle the dataset describing solid phase as either precipitated or sorbed species.
 *
 *
 * TODO:
 * sorption.cc:
 * - Why tests from line 151- 162
 * - use just one switch according to isotherm type
 * - what about time dependent sorption parameters?
 * - line 260: existence of appropriate table should be tested, faster and simpler
 *   even better have method is_precomputed() of the Isotherm class.
 *
 * - consider make sorbed concentration internaly also kg (of substnance) /kg (of rock) and
 *   convert it to mol/kg only on output, should be faster and safer
 *
 * - move all code that computes only values of one isotherm (of any type) to some other class
 *   proposed IsothermFactory
 *
 * - Idea is to have one Isotherm factory object, all necessary fields availabel in Sorption::EqDAta object
 *   and have virtual isotherm_reinit method, called during table initialization and/or in iterative solution.
 *
 *   Then we can have various derived Sorption classes for various purposes. Changing just the isotherm_reinit method we can
 *   change the parameters actually used in isotherms.
 *
 *   Need prototype of dual-porosity reaction class to design precise placement of various involved fields.
 *
 */
#ifndef SORPTION_BASE
#define SORPTION_BASE

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
    EqData();

    Field<3, FieldValue<3>::EnumVector > adsorption_type; ///< Discrete need Selection for initialization.
    Field<3, FieldValue<3>::Scalar > rock_density; ///< Rock matrix density.
    Field<3, FieldValue<3>::Vector > isotherm_mult; ///< Multiplication coefficients (k, omega) for all types of isotherms. Langmuir: c_s = omega * (alpha*c_a)/(1- alpha*c_a), Linear: c_s = k*c_a
    Field<3, FieldValue<3>::Vector > isotherm_other; ///< Langmuir sorption coeficients alpha (in fraction c_s = omega * (alpha*c_a)/(1- alpha*c_a)).
    Field<3, FieldValue<3>::Vector> init_conc_solid; ///< Initial sorbed concentrations. 

    Field<3, FieldValue<3>::Scalar > porosity; ///< Porosity field copied from transport
    
    MultiField<3, FieldValue<3>::Scalar>  conc_solid;    ///< Calculated sorbed concentrations, for output only.

    /// Fields indended for output, i.e. all input fields plus those representing solution.
    FieldSet output_fields;

  };

  /**
   *  Constructor with parameter for initialization of a new declared class member
   */
  SorptionBase(Mesh &init_mesh, Input::Record in_rec, vector<string> &names);
  /**
   * Destructor.
   */
  virtual ~SorptionBase(void);
  /**
   * Prepared to compute sorption inside all of considered elements. 
   * It calls compute_reaction(...) for all the elements controled by concrete processor, when the computation is paralelized.
   */
  virtual void update_solution(void);
  
  /**
   * Sets the output names of substances. 
   * This way we do not overwrite the output of substances in transport
   * e.g.:
   * A -> A             (multifield)-> A_sorbed                (sorption in transport)
   * A -> A_mobile      (multifield)-> A_mobile_sorbed         (sorption in dual porosity - mobile)
   * A -> A_immobile    (multifield)-> A_immobile_sorbed       (sorption in dual porosity - immobile)
   */
  virtual void set_output_names(void);
  
  virtual Input::Type::Selection get_output_selection() = 0;

  /**
   * Initialization routines that are done in constructors of descendants.
   * Method data() which access EqData is pure virtual and cannot be called from the base constructor.
   */
  //void data_initialization(void);
  /**
   * Sets porosity field - makes a field copy from transport.
   */
  void set_porosity(Field<3, FieldValue<3>::Scalar > &por_m);
  
  /**
   * Creates interpolation table for isotherms.
   */
  void make_tables(void);
  
  /// Initializes private members of sorption from the input record.
  void init_from_input(Input::Record in_rec) override;

  void initialize(OutputTime *stream) override;
  void output_data(void) override;
  void output_vector_gather(void) override;
  
  /**
   * Meaningless inherited method.
   */
  void set_concentration_vector(Vec &vec) override;
    
protected:
  /**
   * This method disables to use constructor without parameters.
   */
  SorptionBase();
  
  void initialize_substance_ids(const std::vector<string> &names, Input::Record in_rec);
  
  /** Initializes possible following reactions from input record.
   * It should be called after setting mesh, time_governor, distribution and concentration_matrix
   * if there are some setting methods for reactions called (they are not at the moment, so it could be part of init_from_input).
   */
  void init_from_input_reaction(Input::Record in_rec);
  
  /**
   * For simulation of sorption in just one element either inside of MOBILE or IMMOBILE pores.
   */
  double **compute_reaction(double **concentrations, int loc_el);
  /**
   *
   */
  virtual void isotherm_reinit(std::vector<Isotherm> &isotherms, const ElementAccessor<3> &elm) = 0;
  
  /**
   * or printing parameters of isotherms under consideration, not necessary to store
   */
  void print_sorption_parameters(void);
  
  void allocate_output_mpi(void);
  
  /// Equation field data.
  virtual EqData &data() = 0;

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
  
  unsigned int n_substances_;   //< number of substances that take part in the sorption model
  
  /// Output names of substances and fields respectively.
  std::vector<std::string> output_names_;
  /**
   * Input data set - fields in this set are read from the input file.
   */
  FieldSet input_data_set_;
  /**
   * Array for storage infos about sorbed species concentrations.
   */
  double** conc_solid;
  
  Input::Array input_data;

  Input::Array output_array;

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

#endif
