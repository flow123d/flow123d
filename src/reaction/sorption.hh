/** @brief class Sorption is used to enable simulation of sorption described by either linear or Langmuir isotherm in combination with limited solubility under consideration.
 *
 * Class in this file makes it possible to handle the dataset describing solid phase as either precipitated or sorbed species.
 *
 */
#ifndef SORPTION
#define SORPTION

#include <vector>
#include <input/input_type.hh>

#include "fields/field_base.hh"
#include "reaction/isotherm.hh"

class Mesh;
class Distribution;
class Reaction;
class Isotherm;

/*enum SorptionType {
	none = 0,
	linear = 1,
	langmuir = 2,
	freundlich = 3
};*/

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
			 * 	Sorption type specifies a kind of isothermal description of adsorption.
			 */
			static Input::Type::Selection sorption_type_selection;

			/// Collect all fields
			EqData();

			/**
			 * Overrides EqDataBase::read_bulk_list_item, implements reading of
			 * - init_piezo_head key
			 */
			//RegionSet read_bulk_list_item(Input::Record rec);

			Field<3, FieldValue<3>::EnumVector > sorption_types; // Discrete need Selection for initialization.
			//Field<3, FieldValue<3>::Vector > sorption_types; // Discrete need Selection for initialization.
			Field<3, FieldValue<3>::Scalar > mob_porosity; // Mobile porosity.
			//Field<3, FieldValue<3>::Scalar > immob_porosity; // Immobile porosity.
			Field<3, FieldValue<3>::Scalar > rock_density; // Rock matrix density.
			Field<3, FieldValue<3>::Vector > mult_coefs; // Multiplication coefficients (k, omega) for all types of isotherms. Langmuir: c_s = omega * (alpha*c_a)/(1- alpha*c_a), Linear: c_s = k*c_a
			Field<3, FieldValue<3>::Vector > second_params; // Langmuir sorption coeficients alpha (in fraction c_s = omega * (alpha*c_a)/(1- alpha*c_a)).
		};
        /**
         *  Constructor with parameter for initialization of a new declared class member
         *  TODO: parameter description
         */
		Sorption(Mesh &init_mesh, Input::Record in_rec, vector<string> &names);
		/**
		*	Destructor.
		*/
		~Sorption(void);
        /**
        *
        */
        //void set_mesh_(Mesh *mesh_in);
		/**
		*	For simulation of sorption in just one element either inside of MOBILE or IMMOBILE pores.
		*/
		//virtual
		double **compute_reaction(double **concentrations, int loc_el);
		/**
		*	Prepared to compute sorption inside all of considered elements. It calls compute_reaction(...) for all the elements controled by concrete processor, when the computation is paralelized.
		*/
		//virtual
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
		*
		*/
		void set_sorb_conc_array(double** sorb_conc_array);
		/**
		* This is the way to get bulk parameters from Transport EqData to those in Sorption class, similar to set_sorption_fields in Semchem_interface
		*/
		void set_sorption_fields(Field<3, FieldValue<3>::Scalar> *por_m);
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
		*	Fuctions holds together setting of isotopes, bifurcations and substance indices.
		*/
		void prepare_inputs(Input::Record in_rec);
		/**
		*	For printing parameters of isotherms under consideration, not necessary to store
		*/
		void print_sorption_parameters(void);
		/**
		*
		*/
		//Mesh *mesh_;
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
		* 	Indentifier of the region where sorption take place. region_id
		*/
		std::vector<unsigned int> region_ids;
		/**
		* 	Density of the rock-matrix. Depends on region.
		*/
		//std::vector<double> rock_dens;
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
		* 	Critical concentrations of species dissolved in water.
		*/
		std::vector<double> c_aq_max;
		/**
		*	Three dimensional array contains intersections between isotherms and mass balance lines. It describes behaviour of sorbents in mobile pores of various rock matrix enviroments.
		*	 Up to |nr_of_region x nr_of_substances x n_points| doubles. Because of equidistant step lenght in cocidered system of coordinates, just function values are stored.
		*/
		std::vector<std::vector<Isotherm> > isotherms_mob;
		/**
		*	Three dimensional array contains intersections between isotherms and mass balance lines. It describes behaviour of sorbents in immobile pores of various rock matrix enviroments.
		*	 Up to |nr_of_region x nr_of_substances x n_points| doubles. Because of equidistant step lenght in cocidered system of coordinates, just function values are stored.
		*/
		std::vector<std::vector<Isotherm> > isotherms_immob;
		/**
		* 	Number of points as the base for interpolation.
		*/
		std::vector<std::vector<int> > n_points;
		/**
		* 	Region characteristic inputs.
		*/
		EqData data_;
		/**
		* 	Temporary step_length in rotated system of coordinates. Should be |nr_of_region x nr_of_substances| matrix later.
		*/
		//double step_length;
		/**
		* Array for storage infos about sorbed species concentrations.
		*/
		double** sorbed_conc_array;
	    /**
	     * pointers to sorption fields from transport
	     */
	    Field<3, FieldValue<3>::Scalar > *mob_porosity_;
};

#endif
