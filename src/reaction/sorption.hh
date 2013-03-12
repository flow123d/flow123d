/** @brief class Sorption is used to enable simulation of sorption described by either linear or Langmuir isotherm in combination with limited solubility under consideration.
 *
 * Class in this file makes it possible to handle the dataset describing solid phase as either precipitated or sorbed species.
 *
 */
#ifndef SORPTION
#define SORPTION

#include <vector>
#include <input/input_type.hh>
//#include <reaction/isotherms.hh>

#include "fields/field_base.hh"

class Mesh;
class Distribution;
class Reaction;
class Isotherm;

enum SorptionType {
	none = 0,
	linear = 1,
	langmuir = 2,
	freundlich = 3
};

class Sorption:  public Reaction
{
	public:
		class EqData : public EqDataBase // should be written in class Sorption
		{
		public:
			/**
			 * 	Sorption type specifies a kind of isothermal description of adsorption.
			 */

			static Input::Type::Selection sorption_type_selection;

			/// Collect all fields
			EqData(const std::string &name=0);

			/**
			 * Overrides EqDataBase::read_bulk_list_item, implements reading of
			 * - init_piezo_head key
			 */
			//RegionSet read_bulk_list_item(Input::Record rec);

			Field<3, FieldValue<3>::EnumVector > sorption_types; // Discrete need Selection for initialization.
			Field<3, FieldValue<3>::Scalar > mob_porosity; // Mobile porosity.
			Field<3, FieldValue<3>::Scalar > immob_porosity; // Immobile porosity.
			Field<3, FieldValue<3>::Scalar > rock_density; // Rock matrix density.
			Field<3, FieldValue<3>::Vector > mult_coefs; // Multiplication coefficients (k, omega) for all types of isotherms. Langmuir: c_s = omega * (alpha*c_a)/(1- alpha*c_a), Linear: c_s = k*c_a
			Field<3, FieldValue<3>::Vector > alphas; // Langmuir sorption coeficients alpha (in fraction c_s = omega * (alpha*c_a)/(1- alpha*c_a)).
		};

		/*
	 	* Static variable for new input data types input
		*/
		static Input::Type::Record input_type;
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
		*	For simulation of sorption in just one element either inside of MOBILE or IMMOBILE pores.
		*/
		//virtual
		double **compute_reaction(double **concentrations, int loc_el);
		/**
		*	Prepared to compute sorption inside all of considered elements. It calls compute_reaction(...) for all the elements controled by concrete processor, when the computation is paralelized.
		*/
		//virtual
		virtual void compute_one_step(void);
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
		double step_length;
};

#endif
