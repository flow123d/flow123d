#include <boost/foreach.hpp>

#include "reaction/reaction.hh"
#include "reaction/linear_reaction.hh"
#include "reaction/pade_approximant.hh"
#include "reaction/dual_por_exchange.hh"
#include "semchem/semchem_interface.hh"
#include "reaction/isotherm.hh"
#include "reaction/sorption.hh"

#include "system/system.hh"
#include "system/sys_profiler.hh"

#include "la/distribution.hh"
#include "mesh/mesh.h"
#include "mesh/elements.h"
#include "input/type_selection.hh"

#include "fields/field_set.hh"
#include "fields/field_elementwise.hh" 


using namespace std;
using namespace Input::Type;

Selection SorptionBase::EqData::sorption_type_selection = Selection("AdsorptionType")
	.add_value(Isotherm::none,"none", "No adsorption considered.")
	.add_value(Isotherm::linear, "linear",
			"Linear isotherm runs the concentration exchange between liquid and solid.")
	.add_value(Isotherm::langmuir, "langmuir",
			"Langmuir isotherm runs the concentration exchange between liquid and solid.")
	.add_value(Isotherm::freundlich, "freundlich",
			"Freundlich isotherm runs the concentration exchange between liquid and solid.");



Record SorptionBase::input_type
	= Record("Adsorption", "AUXILIARY RECORD. Should not be directly part of the input tree.")
    .declare_key("substances", Array(String()), Default::obligatory(),
                 "Names of the substances that take part in the adsorption model.")
	.declare_key("solvent_density", Double(), Default("1.0"),
				"Density of the solvent.")
	.declare_key("substeps", Integer(), Default("1000"),
				"Number of equidistant substeps, molar mass and isotherm intersections")
	.declare_key("molar_mass", Array(Double()), Default::obligatory(),
							"Specifies molar masses of all the adsorbing species.")
	.declare_key("solubility", Array(Double(0.0)), Default::optional(), //("-1.0"), //
							"Specifies solubility limits of all the adsorbing species.")
	.declare_key("table_limits", Array(Double(0.0)), Default::optional(), //("-1.0"), //
							"Specifies highest aqueous concentration in interpolation table.")
    .declare_key("input_fields", Array(EqData("").input_data_set_.make_field_descriptor_type("Sorption")), Default::obligatory(), //
                    "Containes region specific data necessary to construct isotherms.")//;
    .declare_key("reaction", ReactionTerm::input_type, Default::optional(), "Reaction model following the sorption.");
    

SorptionBase::EqData::EqData(const string &output_field_name)
{
    ADD_FIELD(rock_density, "Rock matrix density.", "0.0");

    ADD_FIELD(sorption_type,"Considered adsorption is described by selected isotherm."); //
              sorption_type.input_selection(&sorption_type_selection);

    ADD_FIELD(isotherm_mult,"Multiplication parameters (k, omega) in either Langmuir c_s = omega * (alpha*c_a)/(1- alpha*c_a) or in linear c_s = k * c_a isothermal description.","1.0");

    ADD_FIELD(isotherm_other,"Second parameters (alpha, ...) defining isotherm  c_s = omega * (alpha*c_a)/(1- alpha*c_a).","1.0");
    ADD_FIELD(init_conc_solid, "Initial solid concentration of substances."
            " Vector, one value for every substance.", "0");
    
    rock_density.units("");
    init_conc_solid.units("M/L^3");

    input_data_set_ += *this;

    // porosity field is set from governing equation (transport) later
    // hence we do not add it to the input_data_set_
    *this += porosity
             .name("porosity")
             .units("1")
             .just_copy();
    
    output_fields += *this;
    output_fields += conc_solid.name(output_field_name).units("M/L^3");
}


SorptionBase::SorptionBase(Mesh &init_mesh, Input::Record in_rec)//
	: ReactionTerm(init_mesh, in_rec),
	  data_(nullptr)
{
  // creating reaction from input and setting their parameters
  make_reactions(input_record_);
}


SorptionBase::~SorptionBase(void)
{
  if(reaction != nullptr) delete reaction;
  
  if (data_ != nullptr) delete data_;

  VecDestroy(vconc_solid);
  VecDestroy(vconc_solid_out);

  for (unsigned int sbi = 0; sbi < names_.size(); sbi++)
  {
    //no mpi vectors
    xfree(conc_solid[sbi]);
  }
  xfree(conc_solid);
}

void SorptionBase::make_reactions(Input::Record in_rec)
{
  //DBGMSG("SorptionBase init_from_input\n");
  Input::Iterator<Input::AbstractRecord> reactions_it = in_rec.find<Input::AbstractRecord>("reaction");
  if ( reactions_it )
  {
    if (reactions_it->type() == Linear_reaction::input_type ) {
        reaction =  new Linear_reaction(*mesh_, *reactions_it);
    } else
    if (reactions_it->type() == Pade_approximant::input_type) {
        reaction = new Pade_approximant(*mesh_, *reactions_it);
    } else
    if (reactions_it->type() == SorptionBase::input_type ) {
        xprintf(UsrErr, "Sorption model cannot have another descendant sorption model.\n");
    } else
    if (reactions_it->type() == DualPorosity::input_type ) {
        xprintf(UsrErr, "Sorption model cannot have descendant dual porosity model.\n");
    } else
    if (reactions_it->type() == Semchem_interface::input_type )
    {
        xprintf(UsrErr, "Semchem chemistry model is not supported at current time.\n");
    } else
    {
        xprintf(UsrErr, "Unknown reactions type in Sorption model.\n");
    }
  } else
  {
    reaction = nullptr;
  }
}

void SorptionBase::initialize()
{
  //DBGMSG("SorptionBase - initialize.\n");
  ASSERT(distribution != nullptr, "Distribution has not been set yet.\n");
  ASSERT(time_ != nullptr, "Time governor has not been set yet.\n");
  ASSERT(output_stream_,"Null output stream.");
  ASSERT_LESS(0, names_.size());
  
  initialize_substance_ids(names_, input_record_); //computes present substances and sets indices
  init_from_input(input_record_);          //reads non-field data from input
  
  nr_of_regions = mesh_->region_db().bulk_size();
  nr_of_points = input_record_.val<int>("substeps");
  
  //isotherms array resized bellow
  isotherms.resize(nr_of_regions);
  for(int i_reg = 0; i_reg < nr_of_regions; i_reg++)
  {
    isotherms[i_reg].resize(n_substances_);
    for(int i_subst = 0; i_subst < n_substances_; i_subst++)
    {
      isotherms[i_reg][i_subst] = Isotherm();
    }
  }   
  
  //allocating new array for sorbed concentrations
  conc_solid = (double**) xmalloc(names_.size() * sizeof(double*));//new double * [n_substances_];
  conc_solid_out = (double**) xmalloc(names_.size() * sizeof(double*));
  for (unsigned int sbi = 0; sbi < names_.size(); sbi++)
  {
    conc_solid[sbi] = (double*) xmalloc(distribution->lsize() * sizeof(double));//new double[ nr_of_local_elm ];
    conc_solid_out[sbi] = (double*) xmalloc(distribution->size() * sizeof(double));
    //zero initialization of solid concentration for all substances
    for(unsigned int i=0; i < distribution->lsize(); i++)
      conc_solid[sbi][i] = 0;
  }

  allocate_output_mpi();
  
  initialize_fields();
  
  if(reaction != nullptr)
  {
    reaction->names(names_)
      .concentration_matrix(concentration_matrix_, distribution, el_4_loc, row_4_el)
      .set_time_governor(*time_);
  }
}


void SorptionBase::initialize_substance_ids(const vector< string >& names, Input::Record in_rec)
{
  Input::Array substances_array = in_rec.val<Input::Array>("substances");
  unsigned int k, global_idx, i_subst = 0;
  bool found;
  
  for(Input::Iterator<string> spec_iter = substances_array.begin<string>(); spec_iter != substances_array.end(); ++spec_iter, i_subst++)
  {
    //finding the name of a substance in the global array of names
    found = false;
    for(k = 0; k < names.size(); k++)
    {
      if (*spec_iter == names[k]) 
      {
        global_idx = k;
        found = true;
        break;
      }
    }
    
    if(!found)
      xprintf(UsrErr,"Wrong name of %d-th substance - not found in global set of transported substances.\n", i_subst);
    
    //finding the global index of substance in the local array
    found = false;
    for(k = 0; k < substance_global_idx_.size(); k++)
    {
      if(substance_global_idx_[k] == global_idx)
      {
        found = true;
        break;
      }
    }
    
    if(!found)
      substance_global_idx_.push_back(global_idx);

  }  
  n_substances_ = substance_global_idx_.size();
}

void SorptionBase::init_from_input(Input::Record in_rec)
{
  // Common data for all the isotherms loaded bellow
	solvent_density = in_rec.val<double>("solvent_density");

  molar_masses.resize( n_substances_ );
  
	Input::Array molar_mass_array = in_rec.val<Input::Array>("molar_mass");

	if (molar_mass_array.size() == molar_masses.size() )   molar_mass_array.copy_to( molar_masses );
	  else  xprintf(UsrErr,"Number of molar masses %d has to match number of adsorbing species %d.\n", molar_mass_array.size(), molar_masses.size());

	Input::Iterator<Input::Array> solub_iter = in_rec.find<Input::Array>("solubility");
	if( solub_iter )
	{
		solub_iter->copy_to(solubility_vec_);
		if (solubility_vec_.size() != n_substances_)
		{
			xprintf(UsrErr,"Number of given solubility limits %d has to match number of adsorbing species %d.\n", solubility_vec_.size(), n_substances_);
		}
	}else{
		// fill solubility_vec_ with zeros or resize it at least
		solubility_vec_.resize(n_substances_);
	}

	Input::Iterator<Input::Array> interp_table_limits = in_rec.find<Input::Array>("table_limits");
	if( interp_table_limits )
	{
		interp_table_limits->copy_to(table_limit_);
		if (table_limit_.size() != n_substances_)
		{
			xprintf(UsrErr,"Number of given table limits %d has to match number of adsorbing species %d.\n", table_limit_.size(), n_substances_);
		}/**/
	}else{
		// fill table_limit_ with zeros or resize it at least
		table_limit_.resize(n_substances_);
	}
}

void SorptionBase::initialize_fields()
{
  ASSERT(n_substances_ > 0, "Number of substances is wrong, they might have not been set yet.\n");
  data_->sorption_type.n_comp(n_substances_);
  data_->isotherm_mult.n_comp(n_substances_);
  data_->isotherm_other.n_comp(n_substances_);
  data_->init_conc_solid.n_comp(n_substances_);

  // read fields from input file
  data_->input_data_set_.set_input_list(input_record_.val<Input::Array>("input_fields"));

  data_->set_mesh(*mesh_);
  data_->set_limit_side(LimitSide::right);

  //initialization of output
  output_array = input_record_.val<Input::Array>("output_fields");
    //initialization of output
  data_->conc_solid.init(names_);
  data_->conc_solid.set_mesh(*mesh_);
  data_->output_fields.output_type(OutputTime::ELEM_DATA);
  for (int sbi=0; sbi<names_.size(); sbi++)
  {
      // create shared pointer to a FieldElementwise and push this Field to output_field on all regions
      std::shared_ptr<FieldElementwise<3, FieldValue<3>::Scalar> > output_field_ptr(
          new FieldElementwise<3, FieldValue<3>::Scalar>(conc_solid_out[sbi], names_.size(), mesh_->n_elements()));
      data_->conc_solid[sbi].set_field(mesh_->region_db().get_region_set("ALL"), output_field_ptr, 0);
  }
  data_->output_fields.set_limit_side(LimitSide::right);
  output_stream_->add_admissible_field_names(output_array, output_selection);
}


void SorptionBase::zero_time_step()
{
  //DBGMSG("SorptionBase - zero_time_step.\n");
  ASSERT(distribution != nullptr, "Distribution has not been set yet.\n");
  ASSERT(time_ != nullptr, "Time governor has not been set yet.\n");
  ASSERT(output_stream_,"Null output stream.");
  ASSERT_LESS(0, names_.size());
  
  data_->set_time(*time_);
  set_initial_condition();
  make_tables();
    
  // write initial condition
  output_vector_gather();
  data_->output_fields.set_time(*time_);
  data_->output_fields.output(output_stream_);
  
  if(reaction != nullptr)
    reaction->zero_time_step();
}

void SorptionBase::set_initial_condition()
{
  for (unsigned int loc_el = 0; loc_el < distribution->lsize(); loc_el++)
  {
    unsigned int index = el_4_loc[loc_el];
    ElementAccessor<3> ele_acc = mesh_->element_accessor(index);
    arma::vec value = data_->init_conc_solid.value(ele_acc.centre(),
        ele_acc);

    //setting initial solid concentration for substances involved in adsorption
    for (int sbi = 0; sbi < n_substances_; sbi++)
    {
      int subst_id = substance_global_idx_[sbi];
      conc_solid[subst_id][loc_el] = value(sbi);
    }
  }
}


void SorptionBase::update_solution(void)
{
  //DBGMSG("Sorption - update_solution\n");
  data_->set_time(*time_); // set to the last computed time

  // if parameters changed during last time step, reinit isotherms and eventualy 
  // update interpolation tables in the case of constant rock matrix parameters
  if(data_->changed())
    make_tables();
    

  START_TIMER("Sorption");
  for (int loc_el = 0; loc_el < distribution->lsize(); loc_el++)
  {
    compute_reaction(concentration_matrix_, loc_el);
  }
  END_TIMER("Sorption");
  
  if(reaction != nullptr) reaction->update_solution();
}

void SorptionBase::make_tables(void)
{
  ElementAccessor<3> elm;
  BOOST_FOREACH(const Region &reg_iter, this->mesh_->region_db().get_region_set("BULK") )
  {
    int reg_idx = reg_iter.bulk_idx();

    if(data_->is_constant(reg_iter))
    {
      ElementAccessor<3> elm(this->mesh_, reg_iter); // constant element accessor
      isotherm_reinit(isotherms[reg_idx],elm);
      //xprintf(MsgDbg,"parameters are constant\n");
      for(int i_subst = 0; i_subst < n_substances_; i_subst++)
      {
        isotherms[reg_idx][i_subst].make_table(nr_of_points);
      }
    }
  }
}

double **SorptionBase::compute_reaction(double **concentrations, int loc_el) // Sorption simulations are realized just for one element.
{
  //DBGMSG("compute_reaction\n");
    ElementFullIter elem = mesh_->element(el_4_loc[loc_el]);
    double porosity;
    double rock_density;
    Region region = elem->region();
    int reg_id_nr = region.bulk_idx();
    int variabl_int = 0;
    int i_subst, subst_id;

    std::vector<Isotherm> & isotherms_vec = isotherms[reg_id_nr];

    //if(reg_id_nr != 0) cout << "region id is " << reg_id_nr << endl;
    
    // Constant value of rock density and mobile porosity over the whole region => interpolation_table is precomputed
    if (isotherms_vec[0].is_precomputed()) 
    {
      for(i_subst = 0; i_subst < n_substances_; i_subst++)
      {
        subst_id = substance_global_idx_[i_subst];
        //DBGMSG("on s_%d precomputed %d\n",subst_id, isotherms_vec[i_subst].is_precomputed());
      
        Isotherm & isotherm = this->isotherms[reg_id_nr][i_subst];
        isotherm.interpolate((concentration_matrix_[subst_id][loc_el]), conc_solid[subst_id][loc_el]);
      }
    }
    else 
    {
      isotherm_reinit(isotherms_vec, elem->element_accessor());
      
      for(i_subst = 0; i_subst < n_substances_; i_subst++)
      {
        subst_id = substance_global_idx_[i_subst];
        Isotherm & isotherm = this->isotherms[reg_id_nr][i_subst];
        isotherm.compute((concentration_matrix_[subst_id][loc_el]), conc_solid[subst_id][loc_el]);
      }
    }
    
  return concentrations;
}


/**************************************** OUTPUT ***************************************************/

void SorptionBase::allocate_output_mpi(void )
{
    int sbi, n_subst, ierr, rank, np; //, i, j, ph;
    n_subst = names_.size();

    vconc_solid = (Vec*) xmalloc(n_subst * (sizeof(Vec)));
    vconc_solid_out = (Vec*) xmalloc(n_subst * (sizeof(Vec))); // extend to all

    for (sbi = 0; sbi < n_subst; sbi++) {
        ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD,1, distribution->lsize(), mesh_->n_elements(), conc_solid[sbi],
                &vconc_solid[sbi]);
        VecZeroEntries(vconc_solid[sbi]);

        ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1, mesh_->n_elements(), conc_solid_out[sbi], &vconc_solid_out[sbi]);
        VecZeroEntries(vconc_solid_out[sbi]);
    }
}


void SorptionBase::output_vector_gather() 
{
    unsigned int sbi/*, rank, np*/;
    IS is;
    VecScatter vconc_out_scatter;
    //PetscViewer inviewer;

    ISCreateGeneral(PETSC_COMM_SELF, mesh_->n_elements(), row_4_el, PETSC_COPY_VALUES, &is); //WithArray
    VecScatterCreate(vconc_solid[0], is, vconc_solid_out[0], PETSC_NULL, &vconc_out_scatter);
    for (sbi = 0; sbi < names_.size(); sbi++) {
        VecScatterBegin(vconc_out_scatter, vconc_solid[sbi], vconc_solid_out[sbi], INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(vconc_out_scatter, vconc_solid[sbi], vconc_solid_out[sbi], INSERT_VALUES, SCATTER_FORWARD);
    }
    //VecView(transport->vconc[0],PETSC_VIEWER_STDOUT_WORLD);
    //VecView(transport->vconc_out[0],PETSC_VIEWER_STDOUT_WORLD);
    VecScatterDestroy(&(vconc_out_scatter));
    ISDestroy(&(is));
}


void SorptionBase::output_data(void )
{
    //DBGMSG("Sorption output\n");
    output_vector_gather();

    int ierr, rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if (rank == 0)
    {
      // Register fresh output data
      data_->output_fields.set_time(*time_);
      data_->output_fields.output(output_stream_);
    }

    //it can call only linear reaction which has no output at the moment
    //if(reaction) reaction->output_data();
    
    //for synchronization when measuring time by Profiler
    MPI_Barrier(MPI_COMM_WORLD);
}
