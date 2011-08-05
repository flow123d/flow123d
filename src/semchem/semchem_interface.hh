#ifndef interfaceH
	#define interfaceH

#include "mesh/elements.h"
#include "system/par_distribution.hh"

class Distribution;

class Semchem_interface
{
	public:
		/**
		*	Semchem interface is the tool to call a simulation of chemical reactions as a part of transport model.
		*/
		Semchem_interface(double timeStep, Mesh * mesh, int nrOfSpecies, bool dualPorosity); //(int nrOfElements, double ***ConcentrationMatrix, Mesh *mesh);
		/**
		*
		*/
		void compute_reaction(bool porTyp, ElementIter ppelm, int poradi, double ***conc);
		/**
		*
		*/
		void compute_one_step(void);
		/**
		*
		*/
		void set_chemistry_computation(void);
		/**
		*
		*/
		void set_timestep(double new_timestep);
		/**
		*
		*/
		void set_dual_porosity(void);
		/**
		*
		*/
		void set_nr_of_elements(int nrOfElements);
		/**
		*
		*/
		void set_concentration_matrix(double ***ConcentrationsMatrix, Distribution *conc_distr, int *el_4_loc);
		/**
		*
		*/
		void set_el_4_loc(int *el_for_loc);
		/**
		*
		*/
		void set_mesh_(Mesh *mesh);
		/**
		*
		*/
		bool semchem_on;
		/**
		*
		*/
		bool dual_porosity_on;
		/**
		*
		*/
		int nr_of_elements;
		/**
		*
		*/
		double ***concentration_matrix;
	private:
		/**
		*
		*/
		double time_step;
		/**
		*
		*/
		Mesh *mesh_;
		/**
		*
		*/
		Distribution *distribution;
		/**
		*
		*/
		int *el_4_loc;
};
#endif
