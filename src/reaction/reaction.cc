#include <stdlib.h>

#include "reaction/reaction.hh"
#include "system/global_defs.h"
#include "mesh/mesh.h"
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


double **ReactionTerm::compute_reaction(double **concentrations, int loc_el)
{
  ASSERT(0,"double **ReactionTerm::compute_reaction(double **concentrations, int loc_el)" 
           "needs to be re-implemented in ancestors.\n");
  return concentrations;
}

void ReactionTerm::choose_next_time(void)
{
  ASSERT(0,"ReactionTerm does not change TimeGovernor.\n");
}
