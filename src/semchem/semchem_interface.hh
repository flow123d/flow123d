#ifndef interfaceH
	#define interfaceH

#include "../mesh/elements.h"

class Semchem_interface
{
	public:
		Semchem_interface(void);
		//void compute_reaction(bool porTyp, ElementIter ppelm, int poradi, double **conc_mob_arr, double **conc_immob_arr);
		void compute_reaction(bool porTyp, ElementIter ppelm, int poradi, double ***conc);
		void compute_one_step(bool porTyp, ElementIter ppelm, double ***conc);
		void set_chemistry_computation(void);
		void set_timestep(double new_timestep);
		bool semchem_on;
	private:
		double time_step;
//		void priprav(void);
//		double change_time_step;
};
#endif
