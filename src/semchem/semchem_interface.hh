#ifndef interfaceH
	#define interfaceH

#include "../mesh/elements.h"

class Semchem_interface
{
	public:
		Semchem_interface(void);
		void compute_reaction(bool porTyp, double time_step, ElementIter ppelm, int poradi, double **conc_mob_arr, double **conc_immob_arr);
		void compute_one_step(bool porTyp, double time_step, ElementIter ppelm, double ***conc);
		bool semchem_on;
	private:
		double set_timestep(double new_timestep);
//		void priprav(void);
//		double change_time_step;
};
#endif
