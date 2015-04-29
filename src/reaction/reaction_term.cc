
#include "reaction/reaction_term.hh"
#include "system/global_defs.h"

#include "io/output_time.hh"
#include "mesh/mesh.h"

using namespace Input::Type;
        
AbstractRecord & ReactionTerm::get_input_type() {
	static AbstractRecord type = AbstractRecord("ReactionTerm",
			"Equation for reading information about simple chemical reactions.");
	type.close();
	return type;
}

Record ReactionTerm::input_type_output_record
    = Record("ReactionTermOutput", "Output setting for transport equations.")
        .declare_key("output_stream", OutputTime::input_type, Default::obligatory(),
                        "Parameters of output stream.");

ReactionTerm::ReactionTerm(Mesh &init_mesh, Input::Record in_rec)
    : EquationBase(init_mesh, in_rec),
      concentration_matrix_(nullptr),
      el_4_loc_(nullptr),
      row_4_el_(nullptr),
      distribution_(nullptr),
      output_stream_(nullptr)

{
}

ReactionTerm::~ReactionTerm()
{
}




void ReactionTerm::choose_next_time(void)
{
  ASSERT(0,"ReactionTerm does not change TimeGovernor.\n");
}
