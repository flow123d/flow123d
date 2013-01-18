/** @brief class Sorption is used to enable simulation of sorption described by either linear or Langmuir isotherm in combination with limited solubility under consideration.
 *
 * Class in this file makes it possible to handle the dataset describing solid phase as either precipitated or sorbed species.
 *
 */
#ifndef SORPTION
#define SORPTION

#include <vector>
#include <input/input_type.hh>
#include <input/accessors.hh>
#include <reaction/isotherms.hh>

class Mesh;
class Distribution;
class Reaction;

class Sorption: public Reaction
{
	public:
		/*
	 	* Static variable for new input data types input
		*/
		static Input::Type::Record input_type;
		/*
	 	* Static variable gets information about particular sorption parameters in selected region
		*/
		static Input::Type::Record input_type_isotherm;
        /**
         *  Constructor with parameter for initialization of a new declared class member
         *  TODO: parameter description
         */
		Sorption(Mesh &init_mesh, MaterialDatabase &material_database, Input::Record in_rec, vector<string> &names);
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
		/**
		*	Following time_step setting methods are obsolete for computation of equilibrial sorption.
		*/
		//void set_time_step(double new_timestep, Input::Record in_rec);
		//virtual void set_time_step(Input::Record in_rec);
		//virtual
		//void set_time_step(double time_step);
	protected:

		/**
		*	This method disables to use constructor without parameters.
		*/
		Sorption();
		/**
		*	Fuction reads necessery informations to describe sorption and to set substance indices.
		*/
		void prepare_inputs(Input::Record in_rec);
		/**
		*	For printing indices of species which sorbe.
		*/
		void print_indices(int dec_nr, int n_subst);
		/**
		*	For printing parameters of isotherms under consideration.
		*/
		void print_sorption_parameters(int n_subst);
		/**
		*	Function determines intersections between an isotherm and conservation of mass describing lines.
		*/
		void determine_crossections(int k_points);
		/**
		* 	Rotates either intersections or all the [c_a,c_s] points around origin.
		*/
		void rotate_points(double angle, double **points);
		/**
		* 	Makes projection of rotated datapoints on rotated isotherm. Use interpolation.
		*/
		void interpolate_datapoints(void);
		/**
		*	Sequence of integers describing an order of substances participating sorption.
		*/
		std::vector <unsigned int> substance_ids;
		/**
		* 	Critical concentrations solvable in water.
		*/
		std::vector <double> c_aq_max;
		/**
		* 	Area identifiers where sorptions take place
		*/
		std::vector <unsigned int> areas;
		/**
		* 	Types of predefined isotherms.
		*/
		std::vector<string> types;
		/**
		*	Linear isotherm tangential direction, slopes.
		*/
		std::vector<double> directs;
		/**
		* 	Langmuirs' multiplication coefficients and alpha parameters.
		*/
		std::vector<std::vector<double> > coefs;
		/**
		*	Two dimensional array contains intersections between isotherms and mass balance lines.
		*/
		std::vector<std::vector<std::vector<double> > >isotherms;
		
};

#endif
