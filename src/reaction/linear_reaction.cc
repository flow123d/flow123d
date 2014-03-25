#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>

#include "reaction/reaction.hh"
#include "reaction/linear_reaction.hh"
#include "reaction/pade_approximant.hh"
#include "system/system.hh"
#include "system/sys_profiler.hh"

#include "la/distribution.hh"
#include "mesh/mesh.h"

#define MOBILE 0
#define IMMOBILE 1


using namespace Input::Type;

Record Linear_reaction::input_type_one_decay_substep
	= Record("Substep", "Equation for reading information about radioactive decays.")
	.declare_key("parent", String(), Default::obligatory(),
				"Identifier of an isotope.")
    .declare_key("half_life", Double(), Default::optional(),
                "Half life of the parent substance.")
    .declare_key("kinetic", Double(), Default::optional(),
                "Kinetic constants describing first order reactions.")
    .declare_key("products", Array(String()), Default::obligatory(),
				"Identifies isotopes which decays parental atom to.")
	.declare_key("branch_ratios", Array(Double()), Default("1.0"),   // default is one product, with ratio == 1.0
				"Decay chain branching percentage.");


Record Linear_reaction::input_type
	= Record("LinearReactions", "Information for a decision about the way to simulate radioactive decay.")
	.derive_from( Reaction::input_type )
    .declare_key("decays", Array( Linear_reaction::input_type_one_decay_substep ), Default::obligatory(),
                "Description of particular decay chain substeps.");


using namespace std;

Linear_reaction::Linear_reaction(Mesh &init_mesh, Input::Record in_rec, vector<string> &names)//(double timeStep, Mesh * mesh, int nrOfSpecies, bool dualPorosity, Input::Record in_rec) //(double timestep, int nrOfElements, double ***ConvectionMatrix)
      : Reaction(init_mesh, in_rec, names),
      reaction_matrix(nullptr)
{
  prev_conc = new double[ n_all_substances_ ];
	//Input::Array names_array = in_rec.val<Input::Array>("substances");
	//nr_of_species = names_array.size();
	//xprintf(Msg,"nr_of_species is %d\n",nr_of_species);
	//Input::Array dec_array = in_rec.val<Input::Array>("decays");
	//nr_of_decays = dec_array.size();
	//nr_of_isotopes = nr_of_decays + 1;//temporary solution
	//print_half_lives(nr_of_decays);
	//xprintf(Msg,"\n1. Linear_reaction constructor runs.\n");
	//set_indices(in_rec);	//It needs to be called separetelly, earlier.
	//xprintf(Msg,"\n2. Linear_reaction constructor runs.\n");
	//set_half_lives(in_rec);
	//xprintf(Msg,"\n3. Linear_reaction constructor runs.\n");
	//set_bifurcation(in_rec); //this probably fails
	//xprintf(Msg,"\n4. Linear_reaction constructor runs.\n");
  init_from_input(in_rec);
  
  if(!output_rec.is_empty()) xprintf(Warn, "The outpur record in LinearReactions will be ignored.\n");
	//set_time_step(0.5);
}

Linear_reaction::~Linear_reaction()
{
  release_reaction_matrix();
  if(prev_conc != nullptr){
    delete[](prev_conc);
  }
}

void Linear_reaction::initialize(void )
{
  ASSERT(distribution != nullptr, "Distribution has not been set yet.\n");
  ASSERT(time_ != nullptr, "Time governor has not been set yet.\n");
  
  allocate_reaction_matrix();
  modify_reaction_matrix();
}


double **Linear_reaction::allocate_reaction_matrix(void) //reaction matrix initialization
{
	unsigned int rows, cols;

	DBGMSG("We are going to allocate reaction matrix\n");
	if (reaction_matrix == nullptr) reaction_matrix = (double **)xmalloc(n_all_substances_ * sizeof(double*));//allocation section
	for(rows = 0; rows < n_all_substances_; rows++){
		reaction_matrix[rows] = (double *)xmalloc(n_all_substances_ * sizeof(double));
	}
	for(rows = 0; rows < n_all_substances_;rows++){
	 for(cols = 0; cols < n_all_substances_; cols++){
		 if(rows == cols)   reaction_matrix[rows][cols] = 1.0;
		 else           	reaction_matrix[rows][cols] = 0.0;
	 }
	}

	//print_reaction_matrix();
	return reaction_matrix;
}

/*double **Linear_reaction::modify_reaction_matrix(Input::Record in_rec) //prepare the matrix, which describes reactions
{
	int rows,cols, index_par, bif_id;
	int *index_child;
	double rel_step, prev_rel_step;

	half_lives = (double *)xmalloc(nr_of_decays * sizeof(double));
	bifurcation.resize(nr_of_decays);

	if(reaction_matrix == NULL){
		xprintf(Msg,"\nReaction matrix pointer is NULL.\n");
		return NULL;
	}

	set_indices(in_rec);
	set_half_lives(in_rec);
	set_bifurcation(in_rec);

	if(nr_of_decays > 0){
		//pole rozpaduu
		Input::Array decay_array = in_rec.val<Input::Array>("decays");
		int dec_nr = -1;
		//iterator do cyklu
		for(Input::Iterator<Input::Record> dec_it = decay_array.begin<Input::Record>(); dec_it != decay_array.end(); ++dec_it, ++dec_nr)
		{
			index_par = substance_ids[dec_nr][0]; // pole indexu parentu, because indecees in input file run from one whereas indeces in C++ run from ZERO
			//int first_index = substance_ids[0];

			if(cols < (nr_of_isotopes - 1)){
				rel_step = time_step/half_lives[dec_nr];
			}
			reaction_matrix[index_par][index_par] = pow(0.5,rel_step);

			Input::Array prod_array = dec_it->val<Input::Array>("products");

			int i = 0;
			for(Input::Iterator<Input::Array> prod_it = prod_array.begin<Input::Array>(); prod_it != prod_array.end(); ++prod_it, ++i)
			{
				//bif_id = cols;
					reaction_matrix[index_par][substance_ids[dec_nr][i]] += (1 - pow(0.5,rel_step)) * bifurcation[dec_nr][substance_ids[dec_nr][i]];
			}
		}
	}
	print_reaction_matrix();//just for control print
	return reaction_matrix;
}*/

double **Linear_reaction::modify_reaction_matrix(void) //All the parameters are supposed to be known
        {
    unsigned int index_par;
    double rel_step;

    if (reaction_matrix == nullptr) {
        xprintf(Warn, "\nReaction matrix pointer is NULL.\n");
        return nullptr;
    }

    for (unsigned int i_decay = 0; i_decay < half_lives.size(); i_decay++) {
        index_par = substance_ids[i_decay][0];
        rel_step = time_->dt() / half_lives[i_decay];
        reaction_matrix[index_par][index_par] = pow(0.5, rel_step);

        for (unsigned int i_product = 1; i_product < substance_ids[i_decay].size(); ++i_product)
            reaction_matrix[index_par][ substance_ids[i_decay][i_product] ]
                                       = (1 - pow(0.5, rel_step))* bifurcation[i_decay][i_product-1];
    }
    // print_reaction_matrix(); //just for control print
    return reaction_matrix;
}

double **Linear_reaction::compute_reaction(double **concentrations, int loc_el) //multiplication of concentrations array by reaction matrix
{
    unsigned int cols, rows;

    if (reaction_matrix == nullptr) return concentrations;

	for(cols = 0; cols < n_all_substances_; cols++){
		prev_conc[cols] = concentrations[cols][loc_el];
		concentrations[cols][loc_el] = 0.0;
	}

	for(rows = 0; rows < n_all_substances_; rows++){
        for(cols = 0; cols < n_all_substances_; cols++){
            concentrations[rows][loc_el] += prev_conc[cols] * reaction_matrix[cols][rows];
        }
    }

	return concentrations;
}


void Linear_reaction::print_half_lives(int nr_of_substances) {
    int i;

    xprintf(Msg, "\nHalf-lives are defined as:");
    for (i = 0; i < (nr_of_substances - 1); i++) {
        if (i < (nr_of_substances - 2)) //cout << " " << half_lives[i] <<",";
            xprintf(Msg, " %f", half_lives[i]);
        if (i == (nr_of_substances - 2)) //cout << " " << half_lives[i] <<"\n";
            xprintf(Msg, " %f\n", this->half_lives[i]);
    }
}

// TODO: check duplicity of parents
//       raise warning if sum of ratios is not one
void Linear_reaction::init_from_input(Input::Record in_rec)
{
    unsigned int idx;

	Input::Array decay_array = in_rec.val<Input::Array>("decays");

	substance_ids.resize( decay_array.size() );
	half_lives.resize( decay_array.size() );
	bifurcation.resize( decay_array.size() );

	int i_decay=0;
	for (Input::Iterator<Input::Record> dec_it = decay_array.begin<Input::Record>(); dec_it != decay_array.end(); ++dec_it, ++i_decay)
	{
		//half-lives determining part
		Input::Iterator<double> it_hl = dec_it->find<double>("half_life");
		if (it_hl) {
		   half_lives[i_decay] = *it_hl;
		} else {
		   it_hl = dec_it->find<double>("kinetic");
		   if (it_hl) {
			   half_lives[i_decay] = log(2)/(*it_hl);
		   } else {
		    xprintf(UsrErr, "Missing half-life or kinetic in the %d-th reaction.\n", i_decay);
		  }
		}

		//indices determining part
		string parent_name = dec_it->val<string>("parent");
		Input::Array product_array = dec_it->val<Input::Array>("products");
		Input::Array ratio_array = dec_it->val<Input::Array>("branch_ratios"); // has default value [ 1.0 ]

		// substance_ids contains also parent id
		if (product_array.size() > 0)   substance_ids[i_decay].resize( product_array.size()+1 );
		else			xprintf(UsrErr,"Empty array of products in the %d-th reaction.\n", i_decay);


		// set parent index
		idx = find_subst_name(parent_name);
		if (idx < n_all_substances_)	substance_ids[i_decay][0] = idx;
		else                		xprintf(UsrErr,"Wrong name of parent substance in the %d-th reaction.\n", i_decay);

		// set products
		unsigned int i_product = 1;
		for(Input::Iterator<string> product_it = product_array.begin<string>(); product_it != product_array.end(); ++product_it, i_product++)
		{
			idx = find_subst_name(*product_it);
			if (idx < n_all_substances_)   substance_ids[i_decay][i_product] = idx;
			else                    	xprintf(Msg,"Wrong name of %d-th product in the %d-th reaction.\n", i_product-1 , i_decay);
		}

		//bifurcation determining part
        if (ratio_array.size() == product_array.size() )   ratio_array.copy_to( bifurcation[i_decay] );
        else            xprintf(UsrErr,"Number of branches %d has to match number of products %d in the %d-th reaction.\n",
                                       ratio_array.size(), product_array.size(), i_decay);

	}
}

void Linear_reaction::update_solution(void)
{
  DBGMSG("LinearReactions - update solution\n");
  if(time_->is_changed_dt())
  {
    release_reaction_matrix();
    allocate_reaction_matrix();
    modify_reaction_matrix();
  }
    //data_.set_time(*time_); // set to the last computed time
	//if timestep changed then modify_reaction_matrix(), not implemented yet
    //DBGMSG("decay step\n");
    if (reaction_matrix == nullptr)   return;

    START_TIMER("linear reaction step");
	for (unsigned int loc_el = 0; loc_el < distribution->lsize(); loc_el++)
	 {
	 	this->compute_reaction(concentration_matrix, loc_el);
// 	    if (dual_porosity_on == true) {
// 	     this->compute_reaction(concentration_matrix[IMMOBILE], loc_el);
//	    }

	 }
    END_TIMER("linear reaction step");
	 return;
}

void Linear_reaction::release_reaction_matrix(void)
{
	if(reaction_matrix != nullptr)
	{
		for(unsigned int i = 0; i < n_all_substances_; i++)
		{
			if(reaction_matrix[i] != nullptr)
			{
				free(reaction_matrix[i]);
				reaction_matrix[i] = nullptr;
			}
		}
		free(reaction_matrix);
		reaction_matrix = nullptr;
	}
}

void Linear_reaction::print_reaction_matrix(void)
{
	unsigned int cols,rows;

	DBGMSG("r mat: %p\n", reaction_matrix);
	if(reaction_matrix != nullptr){
                if(time_ != NULL)
                  xprintf(Msg,"\ntime_step %f,Reaction matrix looks as follows:\n",time_->dt());
		for(rows = 0; rows < n_all_substances_; rows++){
			for(cols = 0; cols < n_all_substances_; cols++){
					xprintf(Msg,"%f\t",reaction_matrix[rows][cols]);
			}
			xprintf(Msg,"\n");
		}
	}else{
		xprintf(Msg,"\nReaction matrix needs to be allocated.\n");
	}
	return;
}

unsigned int Linear_reaction::find_subst_name(const string &name)
{

    unsigned int k=0;
        for(; k < names_.size(); k++)
                if (name == names_[k]) return k;

        return k;
}
