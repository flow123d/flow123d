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
	/*
	 * Static variable for new input data types input
	 */
	static Input::Type::Record input_type;

	class EqData : public FieldSet // should be written in class Sorption
	{
	public:

		/// Collect all fields
		EqData();

		/// Mass transfer coefficients between mobile and immobile pores.
		Field<3, FieldValue<3>::Vector > alphas;
	};
    	/**
     	* 	Pointer to porosity field from transport
    	*/
    	pScalar porosity_;
	    /**
	    * 	Pointer to porosity field from transport
	    */
	    pScalar immob_porosity_;
	    /**
	    *
        */
		Dual_por_exchange(Mesh &init_mesh, Input::Record in_rec, vector<string> &names);
		/**
		*	Destructor.
		*/
		~Dual_por_exchange(void);
		/**
		*	This method enables to change a data source the program is working with, during simulation.
		*/
		void set_immob_concentration_matrix(double **ConcentrationMatrix, Distribution *conc_distr, int *el_4_loc);
		/**
		*
		*/
		void compute_one_step(void);
		/**
		*
		*/
		//void set_nr_transp(int nr_transp_subst);
		/**
		*
		*/
		void set_porosity(pScalar porosity, pScalar immob_porosity);
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
		* fraction of the mobile porosity and the whole porosity, it was meant to be fraction of the total sorption surface exposed to the mobile zone, in interval (0,1).
		* pointer to phi field from transport
		*/
	    pScalar phi_;
	    /**
	    * mass transfer coefficients between mobile and immobile pores
	    */
	    std::vector<double> alpha_;
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
