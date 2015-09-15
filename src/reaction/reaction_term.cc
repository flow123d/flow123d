
#include "reaction/reaction_term.hh"
#include "system/global_defs.h"

#include "io/output_time.hh"
#include "mesh/mesh.h"

using namespace Input::Type;
        
AbstractRecord & ReactionTerm::get_input_type() {
	return AbstractRecord("ReactionTerm",
			"Equation for reading information about simple chemical reactions.")
			.close();
}

ReactionTerm::ReactionTerm(Mesh &init_mesh, Input::Record in_rec)
    : EquationBase(init_mesh, in_rec),
      concentration_matrix_(nullptr),
      el_4_loc_(nullptr),
      row_4_el_(nullptr),
      distribution_(nullptr)

{
}

ReactionTerm::~ReactionTerm()
{
}




void ReactionTerm::choose_next_time(void)
{
  ASSERT(0,"ReactionTerm does not change TimeGovernor.\n");
}
