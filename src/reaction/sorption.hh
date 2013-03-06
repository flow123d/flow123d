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

			Field<3, FieldValue<3>::EnumVector > sorption_types; // Discrete need Selection for initialization.
			//Field<3, FieldValue<3>::Scalar > nr_of_points; // Number of required mass-balance crossection. away, obsolete
			//Field<3, FieldValue<3>::Scalar > region_ident; // Rock matrix identifier. away, obsolete
			Field<3, FieldValue<3>::Scalar > mob_porosity; // Mobile porosity.
			Field<3, FieldValue<3>::Scalar > immob_porosity; // Immobile porosity.
			Field<3, FieldValue<3>::Scalar > rock_density; // Rock matrix density.
			//Field<3, FieldValue<3>::Vector > specie; // Specie names.
			Field<3, FieldValue<3>::Vector > mult_coefs; // Multiplication coefficients (k, omega) for all types of isotherms. Langmuir: c_s = omega * (alpha*c_a)/(1- alpha*c_a), Linear: c_s = k*c_a
			Field<3, FieldValue<3>::Vector > alphas; // Langmuir sorption coeficients alpha (in fraction c_s = omega * (alpha*c_a)/(1- alpha*c_a)).
		};

		/*
	 	* Static variable for new input data types input
		*/
		static Input::Type::Record input_type;
		/**
		* Static variable for new input data types input, probably obsolete
		*/
		//static Input::Type::Record input_type_isotherm;
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
		* 	Computes coeficients for the matrix mutiplication based rotation of the coordinate system from region affected inputs.
		*/
		void compute_rot_coefs(double porosity, double rock_density, int spec_id);
		/**
		* 	It is used to switch rotation matrix entries for the counterclockwise rotation of the system of coordinates.
		*/
		void switch_rot_coefs(void);
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
		* 	Rotates either intersections or all the [c_a,c_s] points around origin. May be, angle is obsolete
		*/
		std::vector<double> rotate_point(std::vector<double> points);
		/**
		* 	Makes projection of rotated datapoints on rotated isotherm. Uses interpolation.
		*/
		double interpolate_datapoint(std::vector<double> rot_point, int region, int specie);
		/**
		* 	Sets step length  for particular isotherms in rotated coordination system.
		*/
		double set_step_length(void);
		/**
		* 	Sets the entries of rotation matrix. Rotates datapoints. Projects them on isotherm. Rotates them back and scales the result.
		*/
		void handle_datapoints(double rock_density, double porosity, std::vector<double> &prev_conc, std::vector<double> isotherm, int reg_id_nr, int i_subst);
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
		*	Five dimensional array contains intersections between isotherms and mass balance lines. It describes behaviour of sorbents in various rock matrix enviroments.
		*	 Up to |2 (mobile|immobile) x nr_of_region x nr_of_substances x 2 (coordinates) x n_points| doubles.
		*/
		//std::vector<std::vector<std::vector<std::vector<std::vector<double> > > > > isotherm;
		/**
		*	Three dimensional array contains intersections between isotherms and mass balance lines. It describes behaviour of sorbents in mobile pores of various rock matrix enviroments.
		*	 Up to |nr_of_region x nr_of_substances x n_points| doubles. Because of equidistant step lenght in cocidered system of coordinates, just function values are stored.
		*/
		std::vector<std::vector<std::vector<double> > > isotherm;
		/**
		*	Three dimensional array contains intersections between isotherms and mass balance lines. It describes behaviour of sorbents in immobile pores of various rock matrix enviroments.
		*	 Up to |nr_of_region x nr_of_substances x n_points| doubles. Because of equidistant step lenght in cocidered system of coordinates, just function values are stored.
		*	 It is probably obsolete to specify different isotherms for mobile and immobile pores.
		*/
		//std::vector<std::vector<std::vector<double> > > isotherm_immob;
		/**
		* 	Specifies sorption type. Depends on region.
		*/
		//std::vector<std::vector<Sorption_type> > type;
		/**
		* 	Coeficients contained in coordination system rotating matrix.
		*/
		std::vector<double> rot_coefs;
		/**
		* 	Number of points as the base for interpolation. Depends on region.
		*/
		//std::vector<std::vector<int> > n_points;
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
