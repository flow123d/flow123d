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

Record DualPorosity::input_type
        = Record("DualPorosity",
            "Dual porosity model in transport problems.\n"
            "Provides computing the concentration of substances in mobile and immobile zone.\n"
            )
    .derive_from(Reaction::input_type)
    .declare_key("data", Array(DualPorosity::EqData().make_field_descriptor_type("DualPorosity")), Default::obligatory(),
                    "Containes region specific data necessary to construct dual porosity model.")
    
    .declare_key("reactions_mob", Reaction::input_type, Default::optional(), "Reaction model in mobile zone.")
    .declare_key("reactions_immob", Reaction::input_type, Default::optional(), "Reaction model in immobile zone.")
    
    .declare_key("output", Reaction::input_type_output_record.copy_keys(DualPorosity::EqData().output_fields.make_output_field_keys()),
                Default::optional(),
                "Parameters of output stream.");
    
DualPorosity::EqData::EqData()
{
  ADD_FIELD(alpha, "Diffusion coefficient of non-equilibrium linear exchange between mobile and immobile zone (dual porosity)."
            " Vector, one value for every substance.", "0");
  ADD_FIELD(immob_porosity, "Porosity of the immobile zone.", "0");
  ADD_FIELD(init_conc_immobile, "Initial concentration of substances in the immobile zone."
            " Vector, one value for every substance.", "0");
  
  alpha.units("");
  immob_porosity.units("0");
  init_conc_immobile.units("M/L^3");
  
  output_fields += *this;
  output_fields += conc_immobile.name("immobile").units("M/L^3");
}

DualPorosity::DualPorosity(Mesh &init_mesh, Input::Record in_rec, vector<string> &names)
	: Reaction(init_mesh, in_rec, names)
{
    //DBGMSG("DualPorosity - constructor\n");
    
    data_.alpha.n_comp(n_all_substances_);
    data_.init_conc_immobile.n_comp(n_all_substances_);
    
    //setting fields that are set from input file
    input_data_set_+=data_;
    input_data_set_.set_input_list(in_rec.val<Input::Array>("data"));
    
    //creating field for porosity that is set later from the governing equation (transport)
    data_+=(data_.porosity
          .name("porosity")
          .units("0")
         );
    
    data_.set_mesh(init_mesh);
    
    data_.set_limit_side(LimitSide::right);
    
  Input::Iterator<Input::Record> out_rec = in_rec.find<Input::Record>("output");
  //output_rec = in_rec.find<Input::Record>("output");
  if(out_rec) output_rec = *out_rec;
}

DualPorosity::~DualPorosity(void)
{
  if(reaction_mob != nullptr) delete reaction_mob;
  if(reaction_immob != nullptr) delete reaction_immob;
  
  if(!output_rec.is_empty())
  {
    VecDestroy(vconc_immobile);
    VecDestroy(vconc_immobile_out);
  }
  
  for (unsigned int sbi = 0; sbi < n_all_substances_; sbi++) 
  {
      //no mpi vectors
      xfree(conc_immob[sbi]);
  }

  xfree(conc_immob);
}


void DualPorosity::init_from_input(Input::Record in_rec)
{ 
  //DBGMSG("dual_por init_from_input\n");
  
  Input::Iterator<Input::AbstractRecord> reactions_it = in_rec.find<Input::AbstractRecord>("reactions_mob");
  if ( reactions_it ) 
  {
    if (reactions_it->type() == Linear_reaction::input_type ) {
        reaction_mob =  new Linear_reaction(*mesh_, *reactions_it, names_);
                
    } else
    if (reactions_it->type() == Pade_approximant::input_type) {
        reaction_mob = new Pade_approximant(*mesh_, *reactions_it, names_ );
    } else
    if (reactions_it->type() == SorptionBase::input_type ) {
        reaction_mob =  new SorptionMob(*mesh_, *reactions_it, names_);
                
       static_cast<SorptionMob *> (reaction_mob) -> set_porosity(data_.porosity);
       static_cast<SorptionMob *> (reaction_mob) -> set_porosity_immobile(data_.immob_porosity);
                
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
    reaction_mob = nullptr;
  }
  
  reactions_it = in_rec.find<Input::AbstractRecord>("reactions_immob");
  if ( reactions_it ) 
  {
    if (reactions_it->type() == Linear_reaction::input_type ) {
        reaction_immob =  new Linear_reaction(*mesh_, *reactions_it, names_);
                
    } else
    if (reactions_it->type() == Pade_approximant::input_type) {
        reaction_immob = new Pade_approximant(*mesh_, *reactions_it, names_ );
    } else
    if (reactions_it->type() == SorptionBase::input_type ) {
        reaction_immob =  new SorptionImmob(*mesh_, *reactions_it, names_);
        
       static_cast<SorptionImmob *> (reaction_immob) -> set_porosity(data_.porosity);        
       static_cast<SorptionImmob *> (reaction_immob) -> set_porosity_immobile(data_.immob_porosity);
                
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
    reaction_immob = nullptr;
  }
}


void DualPorosity::initialize(void )
{ 
  ASSERT(distribution != nullptr, "Distribution has not been set yet.\n");
  ASSERT(time_ != nullptr, "Time governor has not been set yet.\n");
  
  data_.set_time(*time_);
    
  //allocating memory for immobile concentration matrix
  conc_immob = (double**) xmalloc(n_all_substances_ * sizeof(double*));
  conc_immobile_out = (double**) xmalloc(n_all_substances_ * sizeof(double*));
  for (unsigned int sbi = 0; sbi < n_all_substances_; sbi++)
  {
    conc_immob[sbi] = (double*) xmalloc(distribution->lsize() * sizeof(double));
    conc_immobile_out[sbi] = (double*) xmalloc(distribution->lsize() * sizeof(double));
  }
  //DBGMSG("DualPorosity - init_conc_immobile.\n");
    
  //copied from convection set_initial_condition
  //setting initial condition for immobile concentration matrix
  FOR_ELEMENTS(mesh_, elem)
  {
    if (!distribution->is_local(el_4_loc[elem.index()])) continue;

    unsigned int index = el_4_loc[elem.index()] - distribution->begin();
    ElementAccessor<3> ele_acc = mesh_->element_accessor(elem.index());
    arma::vec value = data_.init_conc_immobile.value(elem->centre(), ele_acc);
        
    for (int sbi=0; sbi < n_all_substances_; sbi++)
    {
      conc_immob[sbi][index] = value(sbi);
    }
  }
  
  if (!output_rec.is_empty())
  {
    //initialization of output
    int rank;
    MPI_Comm_rank(PETSC_COMM_SELF, &rank);
    if (rank == 0)
    {
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
        output_stream = OutputTime::output_stream(output_rec.val<Input::Record>("output_stream"));
    }
    
    allocate_output_mpi();
  
    // write initial condition
    output_vector_gather();
    data_.output_fields.set_time(*time_);
    data_.output_fields.output(output_stream);
  }
  
  // creating reactions from input and setting their parameters
  init_from_input(input_record_);
  
  if(reaction_mob != nullptr)
  { 
    reaction_mob->set_time_governor(*time_);
    reaction_mob->set_concentration_matrix(concentration_matrix, distribution, el_4_loc, row_4_el);
    reaction_mob->initialize();
  }
    
  if(reaction_immob != nullptr) 
  {
    reaction_immob->set_time_governor(*time_);
    reaction_immob->set_concentration_matrix(conc_immob, distribution, el_4_loc, row_4_el);
    reaction_immob->initialize();
  }
}


void DualPorosity::update_solution(void) 
{
  DBGMSG("DualPorosity - update solution\n");
  data_.set_time(*time_);
 
  START_TIMER("dual_por_exchange_step");
  for (unsigned int loc_el = 0; loc_el < distribution->lsize(); loc_el++) 
  {
    compute_reaction(conc_immob, loc_el);
  }
  END_TIMER("dual_por_exchange_step");
  
  if(reaction_mob != nullptr) reaction_mob->update_solution();
  if(reaction_immob != nullptr) reaction_immob->update_solution();
}


double **DualPorosity::compute_reaction(double **concentrations, int loc_el) 
{
  double conc_avg = 0.0;
  unsigned int sbi, sbi_loc;
  double cm, pcm, ci, pci, por_m, por_imm, temp_exp;
   
  ElementFullIter ele = mesh_->element(el_4_loc[loc_el]);
  por_m = data_.porosity.value(ele->centre(),ele->element_accessor());
  por_imm = data_.immob_porosity.value(ele->centre(),ele->element_accessor());
  arma::Col<double> alpha_vec = data_.alpha.value(ele->centre(), ele->element_accessor());
  
  if(time_->dt() >= min_dt)
  {
      //TODO:
      //for (sbi = 0; sbi < n_substances_; sbi++) //over substances involved in dual porosity model
      for (sbi = 0; sbi < n_all_substances_; sbi++) //over all substances
      {
        //sbi_loc = substance_id[sbi];    //mapping to global substance index
                //previous values
                pcm = concentration_matrix[sbi][loc_el];
                pci = conc_immob[sbi][loc_el];

                // ---compute average concentration------------------------------------------
                conc_avg = ((por_m * pcm) + (por_imm * pci)) / (por_m + por_imm);

                if ((conc_avg != 0.0) && (por_imm != 0.0)) {
                        temp_exp = exp(-alpha_vec[sbi] * ((por_m + por_imm) / (por_m * por_imm)) * time_->dt());
                        // ---compute concentration in mobile area-----------------------------------
                        cm = (pcm - conc_avg) * temp_exp + conc_avg;

                        // ---compute concentration in immobile area---------------------------------
                        ci = (pci - conc_avg) * temp_exp + conc_avg;
                        // --------------------------------------------------------------------------
//                         DBGMSG("cm: %f  ci: %f  pcm: %f  pci: %f  conc_avg: %f  alpha: %f  por_m: %f  por_imm: %f  time_dt: %f\n",
//                                 cm, ci, pcm, pci, conc_avg, alpha_vec[sbi], por_m, por_imm, time_->dt());
                        concentration_matrix[sbi][loc_el] = cm;
                        conc_immob[sbi][loc_el] = ci;
                }
        }
  }
  else{
      
      for (sbi = 0; sbi < n_all_substances_; sbi++) {
                //previous values
                pcm = concentration_matrix[sbi][loc_el];
                pci = conc_immob[sbi][loc_el];

                if (por_imm != 0.0) {
                        temp_exp = alpha_vec[sbi]*(pci - pcm);
                        // ---compute concentration in mobile area-----------------------------------
                        cm = temp_exp / por_m + pcm;

                        // ---compute concentration in immobile area---------------------------------
                        ci = -temp_exp / por_imm + pci;
                        // --------------------------------------------------------------------------

                        concentration_matrix[sbi][loc_el] = cm;
                        conc_immob[sbi][loc_el] = ci;
                }
        }
  }
  return conc_immob;
}


void DualPorosity::allocate_output_mpi(void )
{
    int sbi, n_subst, ierr, rank, np; //, i, j, ph;
    n_subst = n_all_substances_;

    MPI_Barrier(PETSC_COMM_WORLD);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &np);

    vconc_immobile = (Vec*) xmalloc(n_subst * (sizeof(Vec)));
    
    // if( rank == 0)
    vconc_immobile_out = (Vec*) xmalloc(n_subst * (sizeof(Vec))); // extend to all


    for (sbi = 0; sbi < n_subst; sbi++) {
        ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD,1, distribution->lsize(), mesh_->n_elements(), conc_immob[sbi],
                &vconc_immobile[sbi]);
        VecZeroEntries(vconc_immobile[sbi]);

        //  if(rank == 0)
        ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1, mesh_->n_elements(), conc_immobile_out[sbi], &vconc_immobile_out[sbi]);
        VecZeroEntries(vconc_immobile_out[sbi]);
    }
}


void DualPorosity::output_vector_gather() 
{
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
  if (!output_rec.is_empty())
  {
    DBGMSG("DualPorosity output\n");
    output_vector_gather();

    // Register fresh output data
    data_.output_fields.set_time(*time_);
    data_.output_fields.output(output_stream);
    //for synchronization when measuring time by Profiler
    MPI_Barrier(MPI_COMM_WORLD);
  }
  
  if(reaction_mob != nullptr) reaction_mob->output_data();
  if(reaction_immob != nullptr) reaction_immob->output_data();
}

