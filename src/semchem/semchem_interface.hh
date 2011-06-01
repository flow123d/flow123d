//---------------------------------------------------------------------------
#ifndef interfaceH
#define interfaceH
//---------------------------------------------------------------------------
#ifndef mobile
	#define mobile 0
#endif
//---------------------------------------------------------------------------
#ifndef immobile
	#define immobile 1
#endif
#include "../mesh/elements.h"

class Semchem_interface
{
	public:
		Semchem_interface(void);
		void compute_reactions(bool porTyp, double time_step, ElementIter ppelm, int poradi, double **conc_mob_arr, double **conc_immob_arr);
		bool semchem_on;
	private:
		double set_timestep(double new_timestep);
//		void priprav(void);
//		double change_time_step;
};

//--------------------------------------------------------------------------
//  GLOBALNI PROMENNE
//---------------------------------------------------------------------------
extern struct TS_prm	G_prm;
extern struct TS_lat 	*P_lat;
extern struct TS_che	*P_che;

#endif
