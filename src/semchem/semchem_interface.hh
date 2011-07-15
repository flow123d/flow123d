#ifndef interfaceH
	#define interfaceH

#include "../mesh/elements.h"

class Semchem_interface
{
	public:
		Semchem_interface(int nrOfElements, double ***ConcentrationMatrix);
		//void compute_reaction(bool porTyp, ElementIter ppelm, int poradi, double **conc_mob_arr, double **conc_immob_arr);
		void compute_reaction(bool porTyp, ElementIter ppelm, int poradi, double ***conc);
		void compute_one_step(void);
		void set_chemistry_computation(void);
		void set_timestep(double new_timestep);
		void set_dual_porosity(void);
		void set_nr_of_elements(int nrOfElements);
		void set_concentration_matrix(double ***ConcentrationsMatrix);
		bool semchem_on;
		bool dual_porosity_on;
		int nr_of_elements;
		double ***concentration_matrix;
	private:
		double time_step;
//		void priprav(void);
//		double change_time_step;
};
#endif
