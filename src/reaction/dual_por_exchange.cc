#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <boost/foreach.hpp>

#include "reaction/reaction.hh"
#include "reaction/dual_por_exchange.hh"
#include "system/system.hh"
#include "system/sys_profiler.hh"
#include <petscmat.h>

#include "la/distribution.hh"
#include "mesh/mesh.h"
#include "mesh/elements.h"
#include "mesh/region.hh"
#include "fields/field_elementwise.hh" 

#include "reaction/sorption_dual.hh"
#include "reaction/sorption_immob.hh"
#include "reaction/sorption_mob.hh"
#include "reaction/linear_reaction.hh"
#include "reaction/pade_approximant.hh"
#include "semchem/semchem_interface.hh"

using namespace Input::Type;
using namespace std;


//it is currently switched of (by "0") until the reference tests are created
const double DualPorosity::min_dt = 0;

Selection DualPorosity::EqData::output_selection
		= EqData().output_fields.make_output_field_selection("DualPorosity_Output")
		.close();

Record DualPorosity::input_type
        = Record("DualPorosity",
            "Dual porosity model in transport problems.\n"
            "Provides computing the concentration of substances in mobile and immobile zone.\n"
            )
    .derive_from(ReactionTerm::input_type)
    .declare_key("input_fields", Array(DualPorosity::EqData().make_field_descriptor_type("DualPorosity")), Default::obligatory(),
                    "Containes region specific data necessary to construct dual porosity model.")
    
    .declare_key("reaction_mobile", ReactionTerm::input_type, Default::optional(), "Reaction model in mobile zone.")
    .declare_key("reaction_immobile", ReactionTerm::input_type, Default::optional(), "Reaction model in immobile zone.")
    
    .declare_key("output_fields", Array(EqData::output_selection),
                Default("conc_immobile"), "List of fields to write to output stream.");
    
DualPorosity::EqData::EqData()
{
  ADD_FIELD(diffusion_rate_immobile, "Diffusion coefficient of non-equilibrium linear exchange between mobile and immobile zone (dual porosity)."
            " Vector, one value for every substance.", "0");
  ADD_FIELD(porosity_immobile, "Porosity of the immobile zone.", "0");
  ADD_FIELD(init_conc_immobile, "Initial concentration of substances in the immobile zone."
            " Vector, one value for every substance.", "0");
  
  diffusion_rate_immobile.units("");
  porosity_immobile.units("0");
  init_conc_immobile.units("M/L^3");
  
  output_fields += *this;
  output_fields += conc_immobile.name("conc_immobile").units("M/L^3");
}

DualPorosity::DualPorosity(Mesh &init_mesh, Input::Record in_rec, vector<string> &names)
	: ReactionTerm(init_mesh, in_rec, names)
{
    //DBGMSG("DualPorosity - constructor\n");
    
    data_.diffusion_rate_immobile.n_comp(n_all_substances_);
    data_.init_conc_immobile.n_comp(n_all_substances_);
    
    //setting fields that are set from input file
    input_data_set_+=data_;
    input_data_set_.set_input_list(in_rec.val<Input::Array>("input_fields"));
    
    //creating field for porosity that is set later from the governing equation (transport)
    data_+=(data_.porosity
          .name("porosity")
          .units("1")
         );

    data_.set_mesh(init_mesh);
    
    data_.set_limit_side(LimitSide::right);

    output_array = in_rec.val<Input::Array>("output_fields");


    //DBGMSG("dual_por init_from_input\n");

    Input::Iterator<Input::AbstractRecord> reactions_it = in_rec.find<Input::AbstractRecord>("reaction_mobile");
    if ( reactions_it )
    {
      if (reactions_it->type() == Linear_reaction::input_type ) {
          reaction_mobile =  new Linear_reaction(*mesh_, *reactions_it, names_);

      } else
      if (reactions_it->type() == Pade_approximant::input_type) {
          reaction_mobile = new Pade_approximant(*mesh_, *reactions_it, names_ );
      } else
      if (reactions_it->type() == SorptionMob::input_type ) {
          reaction_mobile =  new SorptionMob(*mesh_, *reactions_it, names_);
      } else
      if (reactions_it->type() == DualPorosity::input_type ) {
          xprintf(UsrErr, "Dual porosity model cannot have another descendant dual porosity model.\n");
      } else
      if (reactions_it->type() == Semchem_interface::input_type )
      {
          xprintf(UsrErr, "Semchem chemistry model is not supported at current time.\n");
      } else
      {
          xprintf(UsrErr, "Wrong reaction type in DualPorosity model.\n");
      }

    } else
    {
      reaction_mobile = nullptr;
    }

    reactions_it = in_rec.find<Input::AbstractRecord>("reaction_immobile");
    if ( reactions_it )
    {
      if (reactions_it->type() == Linear_reaction::input_type ) {
          reaction_immobile =  new Linear_reaction(*mesh_, *reactions_it, names_);

      } else
      if (reactions_it->type() == Pade_approximant::input_type) {
          reaction_immobile = new Pade_approximant(*mesh_, *reactions_it, names_ );
      } else
      if (reactions_it->type() == SorptionImmob::input_type ) {
          reaction_immobile =  new SorptionImmob(*mesh_, *reactions_it, names_);
      } else
      if (reactions_it->type() == DualPorosity::input_type ) {
          xprintf(UsrErr, "Dual porosity model cannot have another descendant dual porosity model.\n");
      } else
      if (reactions_it->type() == Semchem_interface::input_type )
      {
          xprintf(UsrErr, "Semchem chemistry model is not supported at current time.\n");
      } else
      {
          xprintf(UsrErr, "Unknown reactions type in DualPorosity model.\n");
      }

    } else
    {
      reaction_immobile = nullptr;
    }
}

DualPorosity::~DualPorosity(void)
{
  if(reaction_mobile != nullptr) delete reaction_mobile;
  if(reaction_immobile != nullptr) delete reaction_immobile;

//  if(!output_rec.is_empty())
  {
    VecDestroy(vconc_immobile);
    VecDestroy(vconc_immobile_out);
  }

  for (unsigned int sbi = 0; sbi < n_all_substances_; sbi++)
  {
      //no mpi vectors
      xfree(conc_immobile[sbi]);
  }

  xfree(conc_immobile);
}


void DualPorosity::set_porosity(Field<3, FieldValue<3>::Scalar > &por_m)
{
	data_.set_field(data_.porosity.name(),por_m);

	// assign porosities to reactions
	SorptionMob   *sorption_mob   = dynamic_cast<SorptionMob   *>(reaction_mobile);
	SorptionImmob *sorption_immob = dynamic_cast<SorptionImmob *>(reaction_immobile);

	if (sorption_mob != nullptr)
	{
		sorption_mob->set_porosity(data_.porosity);
		sorption_mob->set_porosity_immobile(data_.porosity_immobile);
	}
	if (sorption_immob != nullptr)
	{
		sorption_immob->set_porosity(data_.porosity);
		sorption_immob->set_porosity_immobile(data_.porosity_immobile);
	}
}


void DualPorosity::initialize(OutputTime *stream)
{
	//DBGMSG("DualPorosity - initialize.\n");
	ASSERT(distribution != nullptr, "Distribution has not been set yet.\n");
	ASSERT(time_ != nullptr, "Time governor has not been set yet.\n");
  
	data_.set_time(*time_);
    
	//allocating memory for immobile concentration matrix
	conc_immobile = (double**) xmalloc(n_all_substances_ * sizeof(double*));
	conc_immobile_out = (double**) xmalloc(n_all_substances_ * sizeof(double*));
	for (unsigned int sbi = 0; sbi < n_all_substances_; sbi++)
	{
		conc_immobile[sbi] = (double*) xmalloc(distribution->lsize() * sizeof(double));
		conc_immobile_out[sbi] = (double*) xmalloc(distribution->size() * sizeof(double));
	}
  
	//DBGMSG("DualPorosity - init_conc_immobile.\n");
	//setting initial condition for immobile concentration matrix
	for (unsigned int loc_el = 0; loc_el < distribution->lsize(); loc_el++)
	{
		unsigned int index = el_4_loc[loc_el];
		ElementAccessor<3> ele_acc = mesh_->element_accessor(index);
		arma::vec value = data_.init_conc_immobile.value(ele_acc.centre(), ele_acc);
        
		for (int sbi=0; sbi < n_all_substances_; sbi++)
		{
			conc_immobile[sbi][loc_el] = value(sbi);
		}
	}
  
    //initialization of output
    allocate_output_mpi();
    output_stream = stream;
	data_.conc_immobile.init(names_);
	data_.conc_immobile.set_mesh(*mesh_);
	data_.output_fields.output_type(OutputTime::ELEM_DATA);

	for (int sbi=0; sbi<n_all_substances_; sbi++)
	{
		// create shared pointer to a FieldElementwise and push this Field to output_field on all regions
		std::shared_ptr<FieldElementwise<3, FieldValue<3>::Scalar> > output_field_ptr(
				new FieldElementwise<3, FieldValue<3>::Scalar>(conc_immobile_out[sbi], n_all_substances_, mesh_->n_elements()));
		data_.conc_immobile[sbi].set_field(mesh_->region_db().get_region_set("ALL"), output_field_ptr, 0);
	}
	data_.output_fields.set_limit_side(LimitSide::right);
	output_stream->add_admissible_field_names(output_array, data_.output_selection);
    
    // write initial condition
    output_vector_gather();
    data_.output_fields.set_time(*time_);
    data_.output_fields.output(output_stream);
  
	// creating reactions from input and setting their parameters
	init_from_input(input_record_);

	if(reaction_mobile != nullptr)
	{
		reaction_mobile->set_time_governor(*time_);
		reaction_mobile->set_concentration_matrix(concentration_matrix, distribution, el_4_loc, row_4_el);
		reaction_mobile->initialize(output_stream);
	}

	if(reaction_immobile != nullptr)
	{
		reaction_immobile->set_time_governor(*time_);
		reaction_immobile->set_concentration_matrix(conc_immobile, distribution, el_4_loc, row_4_el);
		reaction_immobile->initialize(output_stream);
	}
}


void DualPorosity::update_solution(void) 
{
  //DBGMSG("DualPorosity - update solution\n");
  data_.set_time(*time_);
 
  START_TIMER("dual_por_exchange_step");
  for (unsigned int loc_el = 0; loc_el < distribution->lsize(); loc_el++) 
  {
    compute_reaction(conc_immobile, loc_el);
  }
  END_TIMER("dual_por_exchange_step");
  
  if(reaction_mobile != nullptr) reaction_mobile->update_solution();
  if(reaction_immobile != nullptr) reaction_immobile->update_solution();
}


double **DualPorosity::compute_reaction(double **concentrations, int loc_el) 
{
  double conc_avg = 0.0;
  unsigned int sbi, sbi_loc;
  double cm, pcm, ci, pci, por_m, por_imm, temp_exp;
   
  ElementFullIter ele = mesh_->element(el_4_loc[loc_el]);
  por_m = data_.porosity.value(ele->centre(),ele->element_accessor());
  por_imm = data_.porosity_immobile.value(ele->centre(),ele->element_accessor());
  arma::Col<double> diff_vec = data_.diffusion_rate_immobile.value(ele->centre(), ele->element_accessor());
  
  if(time_->dt() >= min_dt)
  {
      //TODO:
      //for (sbi = 0; sbi < n_substances_; sbi++) //over substances involved in dual porosity model
      for (sbi = 0; sbi < n_all_substances_; sbi++) //over all substances
      {
        //sbi_loc = substance_id[sbi];    //mapping to global substance index
                //previous values
                pcm = concentration_matrix[sbi][loc_el];
                pci = conc_immobile[sbi][loc_el];

                // ---compute average concentration------------------------------------------
                conc_avg = ((por_m * pcm) + (por_imm * pci)) / (por_m + por_imm);

                if ((conc_avg != 0.0) && (por_imm != 0.0)) {
                        temp_exp = exp(-diff_vec[sbi] * ((por_m + por_imm) / (por_m * por_imm)) * time_->dt());
                        // ---compute concentration in mobile area-----------------------------------
                        cm = (pcm - conc_avg) * temp_exp + conc_avg;

                        // ---compute concentration in immobile area---------------------------------
                        ci = (pci - conc_avg) * temp_exp + conc_avg;
                        // --------------------------------------------------------------------------
//                         DBGMSG("cm: %f  ci: %f  pcm: %f  pci: %f  conc_avg: %f  diff: %f  por_m: %f  por_imm: %f  time_dt: %f\n",
//                                 cm, ci, pcm, pci, conc_avg, diff_vec[sbi], por_m, por_imm, time_->dt());
                        concentration_matrix[sbi][loc_el] = cm;
                        conc_immobile[sbi][loc_el] = ci;
                }
        }
  }
  else{
      
      for (sbi = 0; sbi < n_all_substances_; sbi++) {
                //previous values
                pcm = concentration_matrix[sbi][loc_el];
                pci = conc_immobile[sbi][loc_el];

                if (por_imm != 0.0) {
                        temp_exp = diff_vec[sbi]*(pci - pcm);
                        // ---compute concentration in mobile area-----------------------------------
                        cm = temp_exp / por_m + pcm;

                        // ---compute concentration in immobile area---------------------------------
                        ci = -temp_exp / por_imm + pci;
                        // --------------------------------------------------------------------------

                        concentration_matrix[sbi][loc_el] = cm;
                        conc_immobile[sbi][loc_el] = ci;
                }
        }
  }
  return conc_immobile;
}


void DualPorosity::allocate_output_mpi(void )
{
  //DBGMSG("DualPorosity - allocate_output_mpi.\n");
    int sbi, n_subst, ierr, rank, np; //, i, j, ph;
    n_subst = n_all_substances_;

    MPI_Barrier(PETSC_COMM_WORLD);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &np);

    vconc_immobile = (Vec*) xmalloc(n_subst * (sizeof(Vec)));
    
    // if( rank == 0)
    vconc_immobile_out = (Vec*) xmalloc(n_subst * (sizeof(Vec))); // extend to all


    for (sbi = 0; sbi < n_subst; sbi++) {
        ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD,1, distribution->lsize(), mesh_->n_elements(), conc_immobile[sbi],
                &vconc_immobile[sbi]);
        VecZeroEntries(vconc_immobile[sbi]);

        //  if(rank == 0)
        ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1, mesh_->n_elements(), conc_immobile_out[sbi], &vconc_immobile_out[sbi]);
        VecZeroEntries(vconc_immobile_out[sbi]);
    }
}


void DualPorosity::output_vector_gather() 
{
  //DBGMSG("DualPorosity - output_vector_gather.\n");
    unsigned int sbi/*, rank, np*/;
    IS is;
    VecScatter vconc_out_scatter;
    //PetscViewer inviewer;

//     MPI_Barrier(PETSC_COMM_WORLD);
//     MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
//     MPI_Comm_size(PETSC_COMM_WORLD, &np);

    
    //ISCreateStride(PETSC_COMM_SELF,mesh_->n_elements(),0,1,&is);
    ISCreateGeneral(PETSC_COMM_SELF, mesh_->n_elements(), row_4_el, PETSC_COPY_VALUES, &is); //WithArray
    VecScatterCreate(vconc_immobile[0], is, vconc_immobile_out[0], PETSC_NULL, &vconc_out_scatter);
    for (sbi = 0; sbi < n_all_substances_; sbi++) {
        VecScatterBegin(vconc_out_scatter, vconc_immobile[sbi], vconc_immobile_out[sbi], INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(vconc_out_scatter, vconc_immobile[sbi], vconc_immobile_out[sbi], INSERT_VALUES, SCATTER_FORWARD);
    }
    //VecView(transport->vconc[0],PETSC_VIEWER_STDOUT_WORLD);
    //VecView(transport->vconc_out[0],PETSC_VIEWER_STDOUT_WORLD);
    VecScatterDestroy(&(vconc_out_scatter));
    ISDestroy(&(is));
}


void DualPorosity::output_data(void )
{
    //DBGMSG("DualPorosity output\n");
    output_vector_gather();

    int ierr, rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if (rank == 0)
    {
      // Register fresh output data
      data_.output_fields.set_time(*time_);
      data_.output_fields.output(output_stream);
    }
    
    //for synchronization when measuring time by Profiler
    MPI_Barrier(MPI_COMM_WORLD);
  
  if(reaction_mobile != nullptr) reaction_mobile->output_data();
  if(reaction_immobile != nullptr) reaction_immobile->output_data();
}

