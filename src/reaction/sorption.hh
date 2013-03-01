/** @brief class Sorption is used to enable simulation of sorption described by either linear or Langmuir isotherm in combination with limited solubility under consideration.
 *
 * Class in this file makes it possible to handle the dataset describing solid phase as either precipitated or sorbed species.
 *
 */
#ifndef SORPTION
#define SORPTION

#include <vector>
#include <input/input_type.hh>
#include <reaction/isotherms.hh>

#include "fields/field_base.hh"

class Mesh;
class Distribution;
class Reaction;
class Isotherm;

enum Sorption_type {
	none = 0,
	Linear = 1,
	Langmuir = 2,
	Freundlich = 3
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

			Field<3, FieldValue<3>::EnumVector > sorption_type; // Discrete need Selection for initialization.
			//Field<3, FieldValue<3>::Scalar > nr_of_points; // Number of required mass-balance crossection.
			//Field<3, FieldValue<3>::Vector > region_ident; // Rock matrix identifier.
			Field<3, FieldValue<3>::Scalar > rock_density; // Rock matrix density.
			Field<3, FieldValue<3>::Vector > scale;  // Scale parameter for linear and Langmuir isotherm
			//Field<3, FieldValue<3>::Scalar > omegas; // Langmuir sorption multiplication parameters.
			Field<3, FieldValue<3>::Vector > alpha; // Langmuir sorption coeficients alpha (in fraction).
			//Field<3, FieldValue<3>::Scalar > solvent_density; // Solvent density.

		};

		/*
	 	* Static variable for new input data types input
		*/
		static Input::Type::Record input_type;
		/**
		* Static variable for new input data types input
		*/
		static Input::Type::Record input_type_isotherm;
		/*
	 	* Static variable gets information about particular sorption parameters in selected region
		*/
		//static Input::Type::Record input_type_isotherm;
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
		 *
		 */
		void precompute_isotherm_tables();
		/**
		*	Fuctions holds together setting of isotopes, bifurcations and substance indices.
		*/
		void prepare_inputs(Input::Record in_rec);
		/**
		* 	Method reads inputs and computes ekvidistant distributed points on all the selected isotherm.
		*/
		void compute_isotherms(Input::Record in_rec);
		/**
		*	For printing parameters of isotherms under consideration, not necessary to store
		*/
		void print_sorption_parameters(void);
		/**
		*	Function determines intersections between an isotherm and conservation of mass describing lines.
		*/
		void determine_crossections(void);
		/**
		* 	Rotates either intersections or all the [c_a,c_s] points around origin.
		*/
		void rotate_point(double angle, std::vector<std::vector<double> > points);
		/**
		* 	Makes projection of rotated datapoints on rotated isotherm. Use interpolation.
		*/
		void interpolate_datapoints(void);
		/**
		* 	Number of regions.
		*/
		int nr_of_regions;
		/**
		* 	Number of substances.
		*/
		int nr_of_substances;
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
		//double solvent_dens;
		/**
		* 	Critical concentrations of species dissolved in water.
		*/
		std::vector<double> c_aq_max;
		/**
		*	Linear isotherm tangential direction, slope. //Up to |nr_of_species x nr_of_regions| parameters. Depends on region.
		*/
		//std::vector<std::vector<double> >directs;
		//double slope;
		/**
		* 	Langmuirs' multiplication coefficient. Depends on region.
		*/
		//double omega;
		/**
		* 	Langmuirs' isotherm alpha parameters. Depends on region.
		*/
		//double alpha;
		/**
		*	Five dimensional array contains intersections between isotherms and mass balance lines. It describes behaviour of sorbents in various rock matrix enviroments. Up to |2 (mobile|immobile) x nr_of_region x nr_of_substances x 2 (coordinates) x n_points| doubles.
		*/
		std::vector<std::vector<std::vector<std::vector<std::vector<double> > > > > isotherm;
		/**
		* 	Specifies sorption type. Depends on region.
		*/
		//std::vector<std::vector<Sorption_type> > type;
		/**
		* 	Information about an angle of rotation of the system of coordinates. |2 (mobile|immobile) x nr_of_region x nr_of_substances| doubles. Radians.
		*/
		std::vector<std::vector<std::vector<double> > > angle;
		/**
		* 	Number of points as the base for interpolation. Depends on region.
		*/
		//std::vector<std::vector<int> > n_points;

		EqData data_;
};

#endif
