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

        
AbstractRecord Reaction::input_type
    = AbstractRecord("Reactions", "Equation for reading information about simple chemical reactions.")
        .declare_key("substances", Array(String()), Default::obligatory(),
                     "Names of the substances that take part in the reaction model.");

Record Reaction::input_type_output_record
    = Record("ReactionOutput", "Output setting for transport equations.")
        .declare_key("output_stream", OutputTime::input_type, Default::obligatory(),
                        "Parameters of output stream.");
//         .declare_key("save_step", Double(0.0), Default::obligatory(),
//                         "Interval between outputs.")
//         .declare_key("output_times", Array(Double(0.0)),
//                         "Explicit array of output times (can be combined with 'save_step'.");

Reaction::Reaction(Mesh &init_mesh, Input::Record in_rec, const  vector<string> &names)
    : EquationBase(init_mesh, in_rec),
      names_(names),
      n_all_substances_ (names.size())
{
  DBGMSG("Reaction constructor.\n");
  initialize_substance_ids(names, in_rec);
  
  // register output vectors
  output_rec = in_rec.val<Input::Record>("output");
}

Reaction::~Reaction()
{
}

void Reaction::initialize_substance_ids(const vector< string >& names, Input::Record in_rec)
{
  Input::Array substances_array = in_rec.val<Input::Array>("substances");
  unsigned int k, idx, i_spec = 0;
  
  for(Input::Iterator<string> spec_iter = substances_array.begin<string>(); spec_iter != substances_array.end(); ++spec_iter, i_spec++)
  {
    //finding name in the global array of names
    for(k = 0; k < names.size(); k++)
    {
      if (*spec_iter == names[k]) 
      {
        idx = k;
        break;
      }
    }
    
    if ((idx < names.size()) && (idx >= 0)) 
    {
      substance_id[i_spec] = idx;       //mapping - if not found, it creates new map
    }
      else    xprintf(UsrErr,"Wrong name of %d-th reaction specie - not found in global set of transported substances.\n", i_spec);
    }
    n_substances_ = substance_id.size();
}

void Reaction::set_concentration_matrix(double **ConcentrationMatrix, Distribution *conc_distr, int *el_4_loc, int *row_4_el)
{
  concentration_matrix = ConcentrationMatrix;
  distribution = conc_distr;
  this->el_4_loc = el_4_loc;
  this->row_4_el = row_4_el;
  return;
}


double **Reaction::compute_reaction(double **concentrations, int loc_el) //multiplication of concentrations array by reaction matrix
{
    cout << "double **Reaction::compute_reaction(double **concentrations, int loc_el) needs to be re-implemented in ancestors." << endl;
        return concentrations;
}

void Reaction::get_parallel_solution_vector(Vec &vec){
	cout << "Reaction.get_parallel_solution_vector(Vec &vec) is not implemented." << endl; //convection->compute_one_step();
}

void Reaction::get_solution_vector(double * &x, unsigned int &a){
	cout << "Reaction.get_solution_vector(double * &x, unsigned int &a) is not implemented." << endl; //convection->compute_one_step();
}

void Reaction::update_solution(void)
{
	cout << "Reaction::update_solution() is not implemented." << endl;
}

void Reaction::choose_next_time(void)
{
	cout << "Reaction::choose_next_time() is not implemented." << endl;
}
