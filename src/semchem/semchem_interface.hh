#ifndef interfaceH
	#define interfaceH

#include "mesh/elements.h"
#include "system/par_distribution.hh"
#include <string.h>

class Distribution;

class Semchem_interface
{
	public:
		/**
		*	Semchem interface is the tool to call a simulation of chemical reactions as a part of transport model. timeStep defines the length of time step for simulation of chemical reactions. nrOfSpecies is the number of transported species. dualPorosity defines type of porosity in examinated soil.
		*/
		Semchem_interface(double timeStep, Mesh * mesh, int nrOfSpecies, bool dualPorosity); //(int nrOfElements, double ***ConcentrationMatrix, Mesh *mesh);
		/**
		*	This method has been prepared to enable simulation of chemical reactions via Semchem. porTyp defines type of porosity. ppelm is a pointer to element we want to simulate chemistry in. poradi is ID of such element. conc is a pointer to threedimensional array full of doubles.
		*/
		void compute_reaction(bool porTyp, ElementIter ppelm, int poradi, double ***conc);
		/**
		*	The function compute_one_step(..) calls compute_reaction(..) for every single element in the mesh.
		*/
		void compute_one_step(void);
		/**
		*	This method reads from ini-file the information if a simulation of chemical raections is switched on.
		*/
		void set_chemistry_computation(void);
		/**
		*	This method enables to change the length of time step for simulation of chemical reactions.
		*/
		void set_timestep(double new_timestep);
		/**
		*	This method reads from ini-file an information if dual porosity is considered in examinated soil.
		*/
		void set_dual_porosity(void);
		/**
		*	This method sets the number of elements contained in mesh.
		*/
		void set_nr_of_elements(int nrOfElements);
		/**
		*	This method sets the pointer to a three dimensional array of doubles.
		*/
		void set_concentration_matrix(double ***ConcentrationsMatrix, Distribution *conc_distr, int *el_4_loc);
		/**
		*	This method sets the pointer to a one dimensional array for converting IDs from local to global.
		*/
		void set_el_4_loc(int *el_for_loc);
		/**
		*	This method sets a pointer to mesh describing examinated area.
		*/
		void set_mesh_(Mesh *mesh);
		/**
		*	function to set  path to an outputfile for semchem-module.
		*/
		void set_fw_chem(std::string semchem_output_file); //(const char* semchem_output_file);
		/**
		*	It containes an information if the simulation of chemical reactions is switched on.
		*/
		bool semchem_on;
		/**
		*	It containes an information if the dual porosity is switched on.
		*/
		bool dual_porosity_on;
		/**
		*	It containes an information about how many elements are contained in mesh.
		*/
		int nr_of_elements;
		/**
		*	It is a pointer on three dimensional matrix full of doubles.
		*/
		double ***concentration_matrix;
		/**
		*	It is name of an output file for semchem.
		*/
		char *fw_chem;
	private:
		/**
		*	It holds an information about the length of time step for chemical reaction simulation.
		*/
		double time_step;
		/**
		*	It is a pointer on mesh.
		*/
		Mesh *mesh_;
		/**
		*	It describes partitioning of elements between processors.
		*/
		Distribution *distribution;
		/**
		*	It enables to change local IDs of elements in mesh into global IDs.
		*/
		int *el_4_loc;
};
#endif
