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

Selection SorptionBase::EqData::sorption_type_selection = Selection("SorptionType")
	.add_value(Isotherm::none,"none", "No adsorption considered")
	.add_value(Isotherm::linear, "linear",
			"Linear isotherm described adsorption considered.")
	.add_value(Isotherm::langmuir, "langmuir",
			"Langmuir isotherm described adsorption considered")
	.add_value(Isotherm::freundlich, "freundlich",
			"Freundlich isotherm described adsorption considered");


Record SorptionBase::input_type
	= Record("Sorption", "Information about all the limited solubility affected adsorptions.")
	.derive_from( Reaction::input_type )
	.declare_key("solvent_dens", Double(), Default("1.0"),
				"Density of the solvent.")
	.declare_key("substeps", Integer(), Default("1000"),
				"Number of equidistant substeps, molar mass and isotherm intersections")
	.declare_key("molar_masses", Array(Double()), Default::obligatory(),
							"Specifies molar masses of all the sorbing species")
	.declare_key("solubility", Array(Double(0.0)), Default::optional(), //("-1.0"), //
							"Specifies solubility limits of all the sorbing species")
	.declare_key("table_limits", Array(Double(0.0)), Default::optional(), //("-1.0"), //
							"Specifies highest aqueous concentration in interpolation table.")
    .declare_key("data", Array(SorptionBase::EqData().make_field_descriptor_type("Sorption")), Default::obligatory(), //
                    "Containes region specific data necessary to construct isotherms.")//;
        
    .declare_key("reactions", Reaction::input_type, Default::optional(), "Reaction model following the sorption.")
    
    .declare_key("output", Reaction::input_type_output_record.copy_keys(SorptionBase::EqData().output_fields.make_output_field_keys()),
                     Default::optional(), "Parameters of output stream.");

SorptionBase::EqData::EqData()
{
    ADD_FIELD(rock_density, "Rock matrix density.", "0.0");

    ADD_FIELD(sorption_types,"Considered adsorption is described by selected isotherm."); //
              sorption_types.input_selection(&sorption_type_selection);

    ADD_FIELD(mult_coefs,"Multiplication parameters (k, omega) in either Langmuir c_s = omega * (alpha*c_a)/(1- alpha*c_a) or in linear c_s = k * c_a isothermal description.","1.0");

    ADD_FIELD(second_params,"Second parameters (alpha, ...) defining isotherm  c_s = omega * (alpha*c_a)/(1- alpha*c_a).","1.0");
    
    rock_density.units("");
    
    output_fields += *this;
    output_fields += conc_sorbed.name("sorbed").units("M/L^3");
}


SorptionBase::SorptionBase(Mesh &init_mesh, Input::Record in_rec, vector<string> &names)//
	: Reaction(init_mesh, in_rec, names)
{
  DBGMSG("SorptionBase constructor.\n");
  
  nr_of_regions = init_mesh.region_db().bulk_size();
  nr_of_points = in_rec.val<int>("substeps");
  
  data_.sorption_types.n_comp(n_substances_);
  data_.mult_coefs.n_comp(n_substances_);
  data_.second_params.n_comp(n_substances_);
    
  data_.set_mesh(init_mesh);
  data_.set_input_list(in_rec.val<Input::Array>("data"));
  
  data_+=(data_.porosity
          .name("porosity")
          .units("0")
         );
  
  data_.set_limit_side(LimitSide::right);
  
  Input::Iterator<Input::Record> out_rec = in_rec.find<Input::Record>("output");
  //output_rec = in_rec.find<Input::Record>("output");
  if(out_rec) output_rec = *out_rec;
  
  output_names_.resize(names_.size());
  
  //Simple vectors holding  common informations.
  molar_masses.resize( n_substances_ );

  //isotherms array resized bellow
  isotherms.resize(nr_of_regions);
  for(int i_reg = 0; i_reg < nr_of_regions; i_reg++)
    for(int i_spec = 0; i_spec < n_substances_; i_spec++)
    {
      Isotherm iso_mob;
      isotherms[i_reg].push_back(iso_mob);
    }
    
  init_from_input(in_rec);
}


SorptionBase::~SorptionBase(void)
{
  if(reaction != nullptr) delete reaction;
  
  if(!output_rec.is_empty())
  {
    VecDestroy(vconc_sorbed);
    VecDestroy(vconc_sorbed_out);
  }

  for (unsigned int sbi = 0; sbi < n_substances_; sbi++) 
  {
    //no mpi vectors
    xfree(sorbed_conc_array[sbi]);
  }
  xfree(sorbed_conc_array);
}


void SorptionBase::init_from_input(Input::Record in_rec)
{ 
    // Common data for all the isotherms loaded bellow
	solvent_dens = in_rec.val<double>("solvent_dens");

	Input::Array molar_mass_array = in_rec.val<Input::Array>("molar_masses");
  
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

void SorptionBase::initialize(void )
{
  ASSERT(distribution != nullptr, "Distribution has not been set yet.\n");
  ASSERT(time_ != nullptr, "Time governor has not been set yet.\n");
  
    //allocating new array for sorbed concentrations
    unsigned int nr_of_local_elm = distribution->lsize();
    sorbed_conc_array = (double**) xmalloc(n_all_substances_ * sizeof(double*));//new double * [n_substances_];
    conc_sorbed_out = (double**) xmalloc(n_all_substances_ * sizeof(double*));
    for (unsigned int sbi = 0; sbi < n_substances_; sbi++)
    {
      sorbed_conc_array[sbi] = (double*) xmalloc(distribution->lsize() * sizeof(double));//new double[ nr_of_local_elm ];
      conc_sorbed_out[sbi] = (double*) xmalloc(distribution->lsize() * sizeof(double));
      for (unsigned int i = 0; i < nr_of_local_elm; i++)
      {
        sorbed_conc_array[sbi][i] = 0.0;
        //conc_sorbed_out[sbi][i] = 0.0;
      }
    }
  
  //initialization of output
  if (!output_rec.is_empty())
  {
    int rank;
    MPI_Comm_rank(PETSC_COMM_SELF, &rank);
    if (rank == 0)
    {
        set_output_names();
        data_.conc_sorbed.init(output_names_);
        data_.conc_sorbed.set_mesh(*mesh_);
        data_.output_fields.output_type(OutputTime::ELEM_DATA);

        for (int sbi=0; sbi<n_all_substances_; sbi++)
        {
                // create shared pointer to a FieldElementwise and push this Field to output_field on all regions
                std::shared_ptr<FieldElementwise<3, FieldValue<3>::Scalar> > output_field_ptr(
                      new FieldElementwise<3, FieldValue<3>::Scalar>(conc_sorbed_out[sbi], n_all_substances_, mesh_->n_elements()));
                data_.conc_sorbed[sbi].set_field(mesh_->region_db().get_region_set("ALL"), output_field_ptr, 0);
        }
        data_.output_fields.set_limit_side(LimitSide::right);
        OutputTime::output_stream(output_rec.val<Input::Record>("output_stream"));
      }
    allocate_output_mpi();
  
  
    DBGMSG("Going to write initial condition.\n");
    // write initial condition
    output_vector_gather();
    data_.output_fields.set_time(*time_);
    data_.output_fields.output(output_rec);
  }
  
  // creating reaction from input and setting their parameters
  init_from_input_reaction(input_record_);
  
  if(reaction != nullptr)
  { 
    reaction->set_time_governor(*time_);
    reaction->set_concentration_matrix(concentration_matrix, distribution, el_4_loc, row_4_el);
    reaction->initialize();
  }
}

void SorptionBase::init_from_input_reaction(Input::Record in_rec)
{
  DBGMSG("SorptionBase init_from_input\n");
  Input::Iterator<Input::AbstractRecord> reactions_it = in_rec.find<Input::AbstractRecord>("reactions");
  if ( reactions_it ) 
  {
    if (reactions_it->type() == Linear_reaction::input_type ) {
        reaction =  new Linear_reaction(*mesh_, *reactions_it, names_);
                
    } else
    if (reactions_it->type() == Pade_approximant::input_type) {
        reaction = new Pade_approximant(*mesh_, *reactions_it, names_ );
        
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

void SorptionBase::update_solution(void)
{
  DBGMSG("SorptionSimple - update_solution\n");
  data_.set_time(*time_); // set to the last computed time
  
  // if parameters changed during last time step, reinit isotherms and eventualy 
  // update interpolation tables in the case of constant rock matrix parameters
  if(data_.changed())
    make_tables();
    

  START_TIMER("Sorption");
  for (int loc_el = 0; loc_el < distribution->lsize(); loc_el++)
  {
    compute_reaction(concentration_matrix, loc_el);
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

    if(data_.is_constant(reg_iter))
    {
      ElementAccessor<3> elm(this->mesh_, reg_iter); // constant element accessor
      isotherm_reinit(isotherms[reg_idx],elm);
      xprintf(MsgDbg,"parameters are constant\n");
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

	std::vector<Isotherm> & isotherms_vec = isotherms[reg_id_nr];

    //if(reg_id_nr != 0) cout << "region id is " << reg_id_nr << endl;
    
    // Constant value of rock density and mobile porosity over the whole region => interpolation_table is precomputed
    if (isotherms_vec[0].is_precomputed()) {
    	for(int i_subst = 0; i_subst < n_substances_; i_subst++)
    	{
    		Isotherm & isotherm = this->isotherms[reg_id_nr][i_subst];
    		int subst_id = substance_id[i_subst];
            isotherm.interpolate((concentration_matrix[subst_id][loc_el]), sorbed_conc_array[i_subst][loc_el]);
    	}
    } else {
		isotherm_reinit(isotherms_vec, elem->element_accessor());
    	for(int i_subst = 0; i_subst < n_substances_; i_subst++)
    	{
            Isotherm & isotherm = this->isotherms[reg_id_nr][i_subst];
            int subst_id = substance_id[i_subst];
            isotherm.compute((concentration_matrix[subst_id][loc_el]), sorbed_conc_array[i_subst][loc_el]);
    	}
    }

	return concentrations;
}

void SorptionBase::set_porosity(Field< 3, FieldValue_< 1, 1, double > >& por_m)
{
  data_.set_field(data_.porosity.name(),por_m); 
}

void SorptionBase::set_output_names(void )
{
  //output names of substances are the same
  //output_names_ = names_;
  for(unsigned int i=0; i < n_all_substances_; i++)
  {
    output_names_[i] = names_[i] + "_mobile";
  }
}


void SorptionBase::print_sorption_parameters(void)
{
  DBGMSG("Not implemented.\n");
  xprintf(Msg, "\nSorption parameters are defined as follows:\n");
}

void SorptionBase::set_concentration_vector(Vec &vc)
{
  DBGMSG("Not implemented.\n");      
}


/**************************************** OUTPUT ***************************************************/

void SorptionBase::allocate_output_mpi(void )
{
    int sbi, n_subst, ierr, rank, np; //, i, j, ph;
    n_subst = n_all_substances_;

    MPI_Barrier(PETSC_COMM_WORLD);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &np);

    vconc_sorbed = (Vec*) xmalloc(n_subst * (sizeof(Vec)));
    
    // if( rank == 0)
    vconc_sorbed_out = (Vec*) xmalloc(n_subst * (sizeof(Vec))); // extend to all


    for (sbi = 0; sbi < n_subst; sbi++) {
        ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD,1, distribution->lsize(), mesh_->n_elements(), sorbed_conc_array[sbi],
                &vconc_sorbed[sbi]);
        VecZeroEntries(vconc_sorbed[sbi]);

        //  if(rank == 0)
        ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1, mesh_->n_elements(), conc_sorbed_out[sbi], &vconc_sorbed_out[sbi]);
        VecZeroEntries(vconc_sorbed_out[sbi]);
    }
}


void SorptionBase::output_vector_gather() 
{
    unsigned int sbi/*, rank, np*/;
    IS is;
    VecScatter vconc_out_scatter;
    //PetscViewer inviewer;

    //  MPI_Barrier(PETSC_COMM_WORLD);
/*    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &np);*/

    
    //ISCreateStride(PETSC_COMM_SELF,mesh_->n_elements(),0,1,&is);
    ISCreateGeneral(PETSC_COMM_SELF, mesh_->n_elements(), row_4_el, PETSC_COPY_VALUES, &is); //WithArray
    VecScatterCreate(vconc_sorbed[0], is, vconc_sorbed_out[0], PETSC_NULL, &vconc_out_scatter);
    for (sbi = 0; sbi < n_all_substances_; sbi++) {
        VecScatterBegin(vconc_out_scatter, vconc_sorbed[sbi], vconc_sorbed_out[sbi], INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(vconc_out_scatter, vconc_sorbed[sbi], vconc_sorbed_out[sbi], INSERT_VALUES, SCATTER_FORWARD);
    }
    //VecView(transport->vconc[0],PETSC_VIEWER_STDOUT_WORLD);
    //VecView(transport->vconc_out[0],PETSC_VIEWER_STDOUT_WORLD);
    VecScatterDestroy(&(vconc_out_scatter));
    ISDestroy(&(is));
}


void SorptionBase::output_data(void )
{
  if (!output_rec.is_empty())
  {
    DBGMSG("Sorption output\n");
    output_vector_gather();

    // Register fresh output data
    data_.output_fields.set_time(*time_);
    data_.output_fields.output(output_rec);

    //for synchronization when measuring time by Profiler
    MPI_Barrier(MPI_COMM_WORLD);
  }
}
