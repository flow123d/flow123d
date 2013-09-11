/** @brief class Sorption_dp is used to enable simulation of sorption described by either linear or Langmuir isotherm in combination with limited solubility under consideration.
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

#include "./reaction/sorption.hh"

class Sorption_dp:  public Sorption
{
	public:
	    /**
	    * 	Pointer to porosity field from transport
	    */
	    pScalar immob_porosity_;
	    /**
	    *
        */
		Sorption_dp(Mesh &init_mesh, Input::Record in_rec, vector<string> &names); //, pScalar mob_porosity, pScalar immob_porosity);
		/**
		*	Destructor.
		*/
		~Sorption_dp(void);
		/**
		*	This method enables to change a data source the program is working with, during simulation.
		*/
		void set_immob_concentration_matrix(double **ConcentrationMatrix, Distribution *conc_distr, int *el_4_loc);
		/**
		*
		*/
		void transport_dual_porosity(void);
		/**
		*
		*/
		void set_scales(double &scale_aqua, double &scale_sorbed, double por_m, double por_imm, double phi, double rock_density, double molar_masses);
		/**
		*
		*/
		void set_nr_transp(int nr_transp_subst);
	protected:
		/**
		*	This method disables to use constructor without parameters.
		*/
		Sorption_dp();
		/**
		*	For printing parameters of isotherms under consideration, not necessary to store
		*/
		void print_sorption_parameters(void);
		/**
		*	Pointer to thwodimensional array[species][elements] containing concentrations either in mobile.
		*/
		double **concentration_matrix;
		/**
		*	Pointer to thwodimensional array[species][elements] containing concentrations either in immobile.
		*/
		double **immob_concentration_matrix;
	    /**
		* fraction of the mobile porosity and the whole porosity, it was meant to be fraction of the total sorption surface exposed to the mobile zone, in interval (0,1).
		* pointer to phi field from transport
		*/
	    pScalar phi_;
	    /**
	    * mass transfer coefficients between mobile and immobile pores
	    */
	    double *alpha_;
		/**
		* 	Number of regions.
		*/
		int nr_of_regions;
		/**
		* 	Number of adsorbing substances.
		*/
		int nr_of_substances;
		/**
		* 	Number of transported substances.
		*/
		int nr_transp_subst_;
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
		* 	Critical concentrations of species dissolved in water.
		*/
		std::vector<double> c_aq_max;
		/**
		*	Three dimensional array contains intersections between isotherms and mass balance lines. It describes behaviour of sorbents in mobile pores of various rock matrix enviroments.
		*	 Up to |nr_of_region x nr_of_substances x n_points| doubles. Because of equidistant step lenght in cocidered system of coordinates, just function values are stored.
		*/
		std::vector<std::vector<Isotherm> > isotherms;
		/**
		* 	Number of points as the base for interpolation.
		*/
		std::vector<std::vector<int> > n_points;
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
