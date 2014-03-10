/** @brief class Dual_por_exchange is used to enable simulation of sorption described by either linear or Langmuir isotherm in combination with limited solubility under consideration.
 *
 * Class in this file makes it possible to handle the dataset describing solid phase as either precipitated or sorbed species.
 *
 */
#ifndef DUAL_POROSITY
#define DUAL_POROSITY

#include <vector>
#include <input/input_type.hh>

#include "fields/field_base.hh"
#include "reaction/isotherm.hh"
#include "./reaction/sorption.hh"

class Mesh;
class Distribution;
class Reaction;

typedef Field<3, FieldValue<3>::Scalar > * pScalar;

#include "./reaction/reaction.hh"

class Dual_por_exchange:  public Reaction
{
public:
  /**
   * Static variable for new input data types input
   */
  static Input::Type::Record input_type;

  class EqData : public EqDataBase // should be written in class Sorption
  {
  public:

    /// Collect all fields
    EqData();

    Field<3, FieldValue<3>::Vector > alpha;            ///< Mass transfer coefficients between mobile and immobile pores.
    Field<3, FieldValue<3>::Scalar > immob_porosity;    ///< Immobile porosity
    
    MultiField<3, FieldValue<3>::Scalar> init_conc_immobile; ///< Initial concentrations in the immobile zone. 

    pScalar porosity; ///< Pointer to mobile porosity
  };

		Dual_por_exchange(Mesh &init_mesh, Input::Record in_rec, vector<string> &names);
		/**
		*	Destructor.
		*/
		~Dual_por_exchange(void);
                
		/**
		*
		*/
		void update_solution(void);
		/**
		*
		*/
		void set_nr_transp(int nr_transp_subst);
		/**
		*
		*/
		void set_porosity(pScalar porosity);
		/// Initialize from input interface.
		virtual void init_from_input(Input::Record in_rec);
	protected:
		/**
		*	This method disables to use constructor without parameters.
		*/
		Dual_por_exchange();
		/**
		*	Pointer to thwodimensional array[species][elements] containing concentrations either in mobile.
		*/
		double **concentration_matrix;
		/**
		*	Pointer to thwodimensional array[species][elements] containing concentrations either in immobile.
		*/
		double **immob_concentration_matrix;


		/**
		* 	Number of adsorbing substances.
		*/
		int nr_of_substances;
		/**
		* 	Number of transported substances.
		*/
		//int nr_transp_subst_;
		/**
		*
		*/
		EqData data_;
		/**
		* Array for storage infos about sorbed species concentrations.
		*/
		double** sorbed_conc_array;
};

#endif
