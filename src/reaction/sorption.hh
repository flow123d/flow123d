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
#ifndef SORPTION
#define SORPTION

#include <vector>
#include <input/input_type.hh>

#include "fields/field_base.hh"
#include "reaction/isotherm.hh"
#include "transport/transport.h"

class Mesh;
class Distribution;
class Reaction;
class Isotherm;

typedef Field<3, FieldValue<3>::Scalar > * pScalar;

class Sorption:  public Reaction
{
	public:
		/*
		 * Static variable for new input data types input
		 */
		static Input::Type::Record input_type;

		class EqData : public EqDataBase // should be written in class Sorption
		{
		public:
			/**
			 * 	Sorption type specifies a kind of equilibrial description of adsorption.
			 */
			static Input::Type::Selection sorption_type_selection;

			/// Collect all fields
			EqData();

			/**
			 * Overrides EqDataBase::read_bulk_list_item, implements reading of
			 * - init_piezo_head key
			 */

			Field<3, FieldValue<3>::EnumVector > sorption_types; // Discrete need Selection for initialization.
			Field<3, FieldValue<3>::Scalar > rock_density; // Rock matrix density.
			Field<3, FieldValue<3>::Vector > mult_coefs; // Multiplication coefficients (k, omega) for all types of isotherms. Langmuir: c_s = omega * (alpha*c_a)/(1- alpha*c_a), Linear: c_s = k*c_a
			Field<3, FieldValue<3>::Vector > second_params; // Langmuir sorption coeficients alpha (in fraction c_s = omega * (alpha*c_a)/(1- alpha*c_a)).
			Field<3, FieldValue<3>::Vector > alphas; // Mass transfer coefficients between mobile and immobile pores.
		};
	    /**
	    * 	Pointer to porosity field from transport
	    */
	    //pScalar mob_porosity_;
	    pScalar porosity_;
	    /**
	    * 	Pointer to porosity field from transport
	    */
	    pScalar immob_porosity_;
        /**
         *  Constructor with parameter for initialization of a new declared class member
         *  TODO: parameter description
         */
		Sorption(Mesh &init_mesh, Input::Record in_rec, vector<string> &names); //, pScalar mob_porosity, pScalar immob_porosity);
		/**
		*	Destructor.
		*/
		~Sorption(void);
		/**
		*	For simulation of sorption in just one element either inside of MOBILE or IMMOBILE pores.
		*/
		double **compute_reaction(double **concentrations, int loc_el);
		/**
		*
		*/
		void isotherm_reinit(std::vector<Isotherm> &isotherms, const ElementAccessor<3> &elm);
		/**
		*	Prepared to compute sorption inside all of considered elements. It calls compute_reaction(...) for all the elements controled by concrete processor, when the computation is paralelized.
		*/
		void compute_one_step(void);
		/**
		*	This method enables to change the timestep for computation of simple chemical reactions. It is obsolete bacause of parent class Reaction.
		*/
		void set_time_step(double new_timestep);
		/**
		* Folowing method enabels the timestep for chemistry to have the value written in ini-file.
		*/
		void set_time_step(Input::Record in_rec);
		/**
		* Inherited init_from_input method extension.
		*/
		void init_from_input(Input::Array bulk_list);
		/**
		*	This method enables to change a data source the program is working with, during simulation.
		*/
		void set_concentration_matrix(double **ConcentrationMatrix, Distribution *conc_distr, int *el_4_loc);
		/**
		*
		*/
		void set_sorb_conc_array(unsigned int nr_of_local_elm);
		/**
		*
		*/
		void set_sorb_conc_array(double** sorb_conc_array);
	    /**
	    *
	    */
	    void set_phi(pScalar phi);
		/**
		* This is the way to get bulk parameters from Transport EqData to those in Sorption_dp class, similar to set_sorption_fields in Semchem_interface
		*/
		void set_porosity(pScalar por, pScalar por_imm);
		/**
		*	Fuctions holds together setting of isotopes, bifurcations and substance indices.
		*/
		void prepare_inputs(Input::Record in_rec, int porosity_type);
		/**
		*	This method enables to change a data source the program is working with, during simulation.
		*/
		void set_immob_concentration_matrix(double **ConcentrationMatrix, Distribution *conc_distr, int *el_4_loc);
		/**
		* Meaningless inherited methods.
		*/
		virtual void update_solution(void);
		virtual void choose_next_time(void);
		virtual void set_time_step_constrain(double dt);
		virtual void get_parallel_solution_vector(Vec &vc);
		virtual void get_solution_vector(double* &vector, unsigned int &size);
	protected:
		/**
		*	This method disables to use constructor without parameters.
		*/
		Sorption();
		/**
		*	For printing parameters of isotherms under consideration, not necessary to store
		*/
		void print_sorption_parameters(void);
		/**
		*	Pointer to thwodimensional array[species][elements] containing concentrations either in mobile or immobile zone.
		*/
		double **concentration_matrix;
		/**
		*	Pointer to thwodimensional array[species][elements] containing concentrations either in immobile.
		*/
		double **immob_concentration_matrix;
	    /**
	    * Mass transfer coefficients between mobile and immobile pores
	    */
	    std::vector<double> alpha_;
		/**
		* 	Number of regions.
		*/
		int nr_of_regions;
		/**
		* 	Number of substances.
		*/
		int nr_of_substances;
		/**
		* 	Temporary nr_of_points can be computed using step_length. Should be |nr_of_region x nr_of_substances| matrix later.
		*/
		int nr_of_points;
		/**
		*	Identifier of the substance undergoing sorption.
		*/
		std::vector<unsigned int> substance_ids;
		/**
		* 	Molar masses of dissolved species (substances)
		*/
		std::vector<double> molar_masses;
		/**
		* 	Density of the solvent. Could be done region dependent, easily.
		*/
		double solvent_dens;
	    /**
		* fraction of the mobile porosity and the whole porosity, it was meant to be fraction of the total sorption surface exposed to the mobile zone, in interval (0,1).
		* pointer to phi field from transport
		*/
	    pScalar phi_;
		/**
		* 	Critical concentrations of species dissolved in water.
		*/
		std::vector<double> c_aq_max;
		/**
		*	Three dimensional array contains intersections between isotherms and mass balance lines. It describes behaviour of sorbents in mobile pores of various rock matrix enviroments.
		*	 Up to |nr_of_region x nr_of_substances x n_points| doubles. Because of equidistant step lenght in cocidered system of coordinates, just function values are stored.
		*/
		std::vector<std::vector<Isotherm> > isotherms;
		/**
		* 	Region characteristic inputs.
		*/
		EqData data_;
		/**
		* Array for storage infos about sorbed species concentrations.
		*/
		double** sorbed_conc_array;
};

#endif
