#include "reaction/reaction_term.hh"
#include "system/global_defs.h"
#include "mesh/mesh.h"
#include "io/output_data.hh"

using namespace Input::Type;
        
AbstractRecord ReactionTerm::input_type
    = AbstractRecord("ReactionTerm", "Equation for reading information about simple chemical reactions.");

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
