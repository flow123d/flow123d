#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>

#include "reaction/reaction.hh"
#include "system/system.hh"
#include "mesh/mesh.h"
#include "mesh/elements.h"
#include "io/output.h"

using namespace Input::Type;
using namespace std;

        
AbstractRecord ReactionTerm::input_type
    = AbstractRecord("ReactionTerm", "Equation for reading information about simple chemical reactions.");

Record ReactionTerm::input_type_output_record
    = Record("ReactionTermOutput", "Output setting for transport equations.")
        .declare_key("output_stream", OutputTime::input_type, Default::obligatory(),
                        "Parameters of output stream.");

ReactionTerm::ReactionTerm(Mesh &init_mesh, Input::Record in_rec)
    : EquationBase(init_mesh, in_rec),
      distribution(nullptr),
      concentration_matrix_(nullptr),
      el_4_loc(nullptr),
      row_4_el(nullptr),
      output_stream_(nullptr)

{
}

ReactionTerm::~ReactionTerm()
{
}



double **ReactionTerm::compute_reaction(double **concentrations, int loc_el) //multiplication of concentrations array by reaction matrix
{
    cout << "double **ReactionTerm::compute_reaction(double **concentrations, int loc_el) needs to be re-implemented in ancestors." << endl;
        return concentrations;
}

void ReactionTerm::get_parallel_solution_vector(Vec &vec){
	cout << "ReactionTerm.get_parallel_solution_vector(Vec &vec) is not implemented." << endl; //convection->compute_one_step();
}

void ReactionTerm::get_solution_vector(double * &x, unsigned int &a){
	cout << "ReactionTerm.get_solution_vector(double * &x, unsigned int &a) is not implemented." << endl; //convection->compute_one_step();
}

void ReactionTerm::choose_next_time(void)
{
	cout << "ReactionTerm::choose_next_time() is not implemented." << endl;
}
