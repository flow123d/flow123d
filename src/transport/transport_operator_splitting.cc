/*
 * transport_operator_splitting.cc
 *
 *  Created on: May 21, 2011
 *      Author: jiri
 */

#include <iostream>
#include <iomanip>

#include "system/system.hh"
#include "system/sys_profiler.hh"
#include "system/xio.h"

#include "transport/transport_operator_splitting.hh"
#include <petscmat.h>
#include "system/sys_vector.hh"
#include "coupling/time_governor.hh"
#include "coupling/equation.hh"
#include "transport/transport.h"
#include "transport/transport_dg.hh"
#include "mesh/mesh.h"
#include "flow/old_bcd.hh"

#include "reaction/reaction.hh"
#include "reaction/linear_reaction.hh"
#include "reaction/pade_approximant.hh"
//#include "reaction/isotherm.hh"
#include "reaction/sorption_dp.hh"

#include "semchem/semchem_interface.hh"

#include "la/distribution.hh"
#include "io/output.h"

#include "input/input_type.hh"
#include "input/accessors.hh"

using namespace Input::Type;

AbstractRecord TransportBase::input_type
	= AbstractRecord("Transport", "Secondary equation for transport of substances.")
	.declare_key("time", TimeGovernor::input_type, Default::obligatory(),
			"Time governor setting for the transport model.")
	.declare_key("substances", Array(String()), Default::obligatory(),
			"Names of transported substances.")
	    // input data
	.declare_key("sorption_enable", Bool(), Default("false"),
			"Model of sorption.")
	.declare_key("dual_porosity", Bool(), Default("false"),
			"Dual porosity model.")
	.declare_key("output", TransportBase::input_type_output_record, Default::obligatory(),
			"Parameters of output stream.");


Record TransportBase::input_type_output_record
	= Record("TransportOutput", "Output setting for transport equations.")
	.declare_key("output_stream", OutputTime::input_type, Default::obligatory(),
			"Parameters of output stream.")
	.declare_key("save_step", Double(0.0), Default::obligatory(),
			"Interval between outputs.")
	.declare_key("output_times", Array(Double(0.0)),
			"Explicit array of output times (can be combined with 'save_step'.")
	.declare_key("conc_mobile_p0", String(),
			"Name of output stream for P0 approximation of the concentration in mobile phase.")
	.declare_key("conc_immobile_p0", String(),
			"Name of output stream for P0 approximation of the concentration in immobile phase.")
	.declare_key("conc_mobile_sorbed_p0", String(),
			"Name of output stream for P0 approximation of the surface concentration of sorbed mobile phase.")
	.declare_key("conc_immobile_sorbed_p0", String(),
			"Name of output stream for P0 approximation of the surface concentration of sorbed immobile phase.");


Record TransportOperatorSplitting::input_type
	= Record("TransportOperatorSplitting",
            "Explicit FVM transport (no diffusion)\n"
            "coupled with reaction and sorption model (ODE per element)\n"
            " via operator splitting.")
    .derive_from(TransportBase::input_type)
	.declare_key("reactions", Reaction::input_type, Default::optional(),
                "Initialization of per element reactions.")
    .declare_key("adsorptions", Sorption::input_type, Default::optional(),
    			"Initialization of per element sorptions.")
    .declare_key("bc_data", Array(ConvectionTransport::EqData().boundary_input_type()
    		.declare_key("old_boundary_file", IT::FileName::input(), "Input file with boundary conditions (obsolete).")
    		.declare_key("bc_times", Array(Double()), Default::optional(),
    				"Times for changing the boundary conditions (obsolete).")
    		), IT::Default::obligatory(), "")
    .declare_key("bulk_data", Array(ConvectionTransport::EqData().bulk_input_type()),
    		IT::Default::obligatory(), "");


TransportBase::TransportEqData::TransportEqData(const std::string& eq_name)
: EqDataBase(eq_name)
{

	ADD_FIELD(init_conc, "Initial concentrations.", Default("0"));
	ADD_FIELD(bc_conc, "Boundary conditions for concentrations.", Default("0"));
	ADD_FIELD(por_m, "Mobile porosity", Default("1"));

	ADD_FIELD(sources_density, "Density of concentration sources.", Default("0"));
	ADD_FIELD(sources_sigma, "Concentration flux.", Default("0"));
	ADD_FIELD(sources_conc, "Concentration sources threshold.", Default("0"));

}


TransportBase::TransportBase(Mesh &mesh, const Input::Record in_rec)
: EquationBase(mesh, in_rec ),
  mh_dh(NULL)
{
	balance_output_file = xfopen( FilePath("mass_balance.txt", FilePath::output_file), "wt");
}

TransportBase::~TransportBase()
{
	if (balance_output_file != NULL) xfclose(balance_output_file);
}


void TransportBase::mass_balance() {
    F_ENTRY;

    // First we compute the quantities, then on process 0 we write the output
    const int nsubst = n_substances(),
    		nboundary = mesh_->region_db().boundary_size(),
    		nbulk = mesh_->region_db().bulk_size();
    vector<vector<double> > bcd_balance(nsubst, vector<double>(nboundary, 0)),
    		bcd_plus_balance(nsubst, vector<double>(nboundary, 0)),
    		bcd_minus_balance(nsubst, vector<double>(nboundary, 0)),
    		mass(nsubst, vector<double>(nbulk, 0)),
    		src_balance(nsubst, vector<double>(nbulk, 0));


    //computing mass fluxes over boundaries and elements
    calc_fluxes(bcd_balance, bcd_plus_balance, bcd_minus_balance);
	calc_elem_sources(mass, src_balance);


	// gather results from processes and sum them up
	int bcd_size = mesh_->region_db().boundary_size();
	int blk_size = mesh_->region_db().bulk_size();
	int buf_size = nsubst*(3*bcd_size + 2*blk_size);
	double sendbuffer[buf_size], recvbuffer[buf_size];

	// prepare sendbuffer
	for (int i=0; i<nsubst; i++)
	{
		for (int j=0; j<bcd_size; j++)
		{
			sendbuffer[i*3*bcd_size+           j] = bcd_balance[i][j];
			sendbuffer[i*3*bcd_size+  bcd_size+j] = bcd_plus_balance[i][j];
			sendbuffer[i*3*bcd_size+2*bcd_size+j] = bcd_minus_balance[i][j];
		}
		for (int j=0; j<blk_size; j++)
		{
			sendbuffer[nsubst*3*bcd_size+i*2*blk_size         +j] = mass[i][j];
			sendbuffer[nsubst*3*bcd_size+i*2*blk_size+blk_size+j] = src_balance[i][j];
		}
	}
	MPI_Reduce(&sendbuffer,recvbuffer,buf_size,MPI_DOUBLE,MPI_SUM,0,PETSC_COMM_WORLD);

	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	if (rank != 0) return;



	// update balance vectors
	for (int i=0; i<nsubst; i++)
	{
		for (int j=0; j<bcd_size; j++)
		{
			bcd_balance[i][j]       = recvbuffer[i*3*bcd_size+           j];
			bcd_plus_balance[i][j]  = recvbuffer[i*3*bcd_size+  bcd_size+j];
			bcd_minus_balance[i][j] = recvbuffer[i*3*bcd_size+2*bcd_size+j];
		}
		for (int j=0; j<blk_size; j++)
		{
			mass[i][j]        = recvbuffer[nsubst*3*bcd_size+i*2*blk_size         +j];
			src_balance[i][j] = recvbuffer[nsubst*3*bcd_size+i*2*blk_size+blk_size+j];
		}
	}


	vector<double> bcd_total_balance(nsubst, 0.), // for computing total balance on boundary
			bcd_total_inflow(nsubst, 0.),
			bcd_total_outflow(nsubst, 0.),
			mass_total(nsubst, 0.),
			src_total_balance(nsubst, 0.);

    //printing the head of mass balance file
    unsigned int c = 6; //column number without label
    unsigned int w = 14;  //column width
    unsigned int wl = 2*(w-5)+7;  //label column width
    stringstream s; //helpful stringstream
    string bc_head_format = "# %-*s%-*s%-*s%-*s%-*s%-*s\n",
           bc_format = "%*s%-*d%-*s%-*s%-*g%-*g%-*g\n",
           bc_total_format = "# %-*s%-*s%-*g%-*g%-*g\n";
    s << setw((w*c+wl-14)/2) << setfill('-') << "--"; //drawing half line
    fprintf(balance_output_file,"# %s MASS BALANCE %s\n",s.str().c_str(), s.str().c_str());
    fprintf(balance_output_file,"# Time: %f\n\n\n",time().t());

    //BOUNDARY
    fprintf(balance_output_file,"# Mass flux through boundary:\n");
    fprintf(balance_output_file,bc_head_format.c_str(),w,"[boundary_id]",wl,"[label]",
                            w,"[substance]",w,"[total flux]",w,"[outward flux]",w,"[inward flux]");
    s.clear();
    s.str(std::string());
    s << setw(w*c+wl) << setfill('-') << "-";
    fprintf(balance_output_file,"# %s\n",s.str().c_str());  //drawing long line


    //printing mass fluxes over boundaries
    DBGMSG("DB[boundary] size: %u\n", mesh_->region_db().boundary_size());
    const RegionSet & b_set = mesh_->region_db().get_region_set("BOUNDARY");

    for( RegionSet::const_iterator reg = b_set.begin(); reg != b_set.end(); ++reg) {
        //DBGMSG("writing reg->idx() and id() and boundary_idx(): %d\t%d\t%d\n", reg->idx(), reg->id(), reg->boundary_idx());
    	for (int sbi=0; sbi<nsubst; sbi++) {
    		bcd_total_balance[sbi] += bcd_balance[sbi][reg->boundary_idx()];
    		bcd_total_outflow[sbi] += bcd_plus_balance[sbi][reg->boundary_idx()];
    		bcd_total_inflow[sbi] += bcd_minus_balance[sbi][reg->boundary_idx()];
			fprintf(balance_output_file, bc_format.c_str(),2,"",w,reg->id(),wl,reg->label().c_str(),
					w, substance_names()[sbi].c_str(),w, bcd_balance[sbi][reg->boundary_idx()],
					w, bcd_plus_balance[sbi][reg->boundary_idx()],
					w, bcd_minus_balance[sbi][reg->boundary_idx()]);
    	}
    }
    //total boundary balance
    fprintf(balance_output_file,"# %s\n",s.str().c_str());  // drawing long line
    for (int sbi=0; sbi<nsubst; sbi++)
    	fprintf(balance_output_file, bc_total_format.c_str(),w+wl,"Total mass flux of substance",
                w,substance_names()[sbi].c_str(),w,bcd_total_balance[sbi], w, bcd_total_outflow[sbi], w, bcd_total_inflow[sbi]);
    fprintf(balance_output_file, "\n\n");


    //SOURCES
    string src_head_format = "# %-*s%-*s%-*s%-*s%-*s\n",
           src_format = "%*s%-*d%-*s%-*s%-*g%-*g\n",
           src_total_format = "# %-*s%-*s%-*g%-*g\n";
    //computing water balance of sources
    fprintf(balance_output_file,"# Mass and sources on regions:\n");   //head
    fprintf(balance_output_file,src_head_format.c_str(),w,"[region_id]",wl,"[label]",
                            w,"[substance]",w,"[total_mass]",w,"[total_source]");
    fprintf(balance_output_file,"# %s\n",s.str().c_str());  //long line

    //printing water balance of sources
    DBGMSG("DB[bulk] size: %u\n", mesh_->region_db().bulk_size());
    const RegionSet & bulk_set = mesh_->region_db().get_region_set("BULK");
    for( RegionSet::const_iterator reg = bulk_set.begin(); reg != bulk_set.end(); ++reg)
    {
    	for (int sbi=0; sbi<nsubst; sbi++)
    	{
    		mass_total[sbi] += mass[sbi][reg->bulk_idx()];
			src_total_balance[sbi] += src_balance[sbi][reg->bulk_idx()];
			//"%*s%-*d%-*s  %-*g%-*s%-*g\n";
			fprintf(balance_output_file, src_format.c_str(), 2,"", w, reg->id(), wl,
					reg->label().c_str(), w, substance_names()[sbi].c_str(),
					w,mass[sbi][reg->bulk_idx()],
					w,src_balance[sbi][reg->bulk_idx()]);
    	}
    }
    // total sources balance
    fprintf(balance_output_file,"# %s\n",s.str().c_str());  //drawing long line
    for (int sbi=0; sbi<nsubst; sbi++)
    	fprintf(balance_output_file, src_total_format.c_str(),w+wl,"Total mass and sources balance",
                w,substance_names()[sbi].c_str(),
                w,mass_total[sbi],
                w,src_total_balance[sbi]);


    // cummulative balance over time

	// quantities need for the balance over time interval
	static vector<double> initial_mass(nsubst, 0.),
			last_sources(nsubst, 0.),
			last_fluxes(nsubst, 0.),
			integrated_sources(nsubst, 0.),
			integrated_fluxes(nsubst, 0.);
	static double initial_time, last_time;
	static bool initial = true;

	if (initial)
	{
		initial_time = time().t();
		last_time = initial_time;
		for (int i=0; i<nsubst; i++)
			initial_mass[i] = mass_total[i];
		initial = false;
	}

	fprintf(balance_output_file, "\n\n# Cumulative mass balance on time interval [%-g,%-g]\n"
			"# %-*s%-*s%-*s%-*s%-*s%-*s%-*s%-*s\n",
			initial_time, time().t(),
			w,"[substance]",
			w,"[A=init. mass]",
			w,"[B=source]",
			w,"[C=flux]",
			w,"[A+B-C]",
			w,"[D=curr. mass]",
			w,"[A+B-C-D=err.]",
			w,"[rel. error]");

	switch (time_scheme())
	{
	case explicit_euler:
		for (int i=0; i<nsubst; i++)
		{
			integrated_sources[i] += last_sources[i]*(time().t()-last_time);
			integrated_fluxes[i] += last_fluxes[i]*(time().t()-last_time);
		}
		break;
	case implicit_euler:
		for (int i=0; i<nsubst; i++)
		{
			integrated_sources[i] += src_total_balance[i]*(time().t()-last_time);
			integrated_fluxes[i] += bcd_total_balance[i]*(time().t()-last_time);
		}
		break;
	case crank_nicholson:
		for (int i=0; i<nsubst; i++)
		{
			integrated_sources[i] += (last_sources[i]+src_total_balance[i])*0.5*(time().t()-last_time);
			integrated_fluxes[i] += (last_fluxes[i]+bcd_total_balance[i])*0.5*(time().t()-last_time);
		}
		break;
	default:
		break;
	}

	for (int i=0; i<nsubst; i++)
	{
		double denominator = max(fabs(initial_mass[i]+integrated_sources[i]-integrated_fluxes[i]),fabs(mass_total[i]));
		fprintf(balance_output_file, "  %-*s%-*g%-*g%-*g%-*g%-*g%-*g%-*g\n",
				w,substance_names()[i].c_str(),
				w,initial_mass[i],
				w,integrated_sources[i],
				w,integrated_fluxes[i],
				w,initial_mass[i]+integrated_sources[i]-integrated_fluxes[i],
				w,mass_total[i],
				w,initial_mass[i]+integrated_sources[i]-integrated_fluxes[i]-mass_total[i],
				w,fabs(initial_mass[i]+integrated_sources[i]-integrated_fluxes[i]-mass_total[i])/(denominator==0?1:denominator));
	}

	last_time = time().t();
	for (int i=0; i<nsubst; i++)
	{
		last_sources[i] = src_total_balance[i];
		last_fluxes[i] = bcd_total_balance[i];
	}

    fprintf(balance_output_file, "\n\n");
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



TransportOperatorSplitting::TransportOperatorSplitting(Mesh &init_mesh, const Input::Record &in_rec)
: TransportBase(init_mesh, in_rec)
{
	time_scheme_ = none;
	Distribution *el_distribution;
	int *el_4_loc;

    // double problem_save_step = OptGetDbl("Global", "Save_step", "1.0");

    in_rec.val<Input::Array>("substances").copy_to(subst_names_);
    n_subst_ = subst_names_.size();
	convection = new ConvectionTransport(*mesh_, in_rec);

	Input::Iterator<Input::AbstractRecord> reactions_it = in_rec.find<Input::AbstractRecord>("reactions");
	if ( reactions_it ) {
		if (reactions_it->type() == Linear_reaction::input_type ) {
	        decayRad =  new Linear_reaction(init_mesh, *reactions_it, subst_names_);
	        convection->get_par_info(el_4_loc, el_distribution);
	        decayRad->set_dual_porosity(convection->get_dual_porosity());
	        static_cast<Linear_reaction *> (decayRad) -> modify_reaction_matrix();
	        decayRad->set_concentration_matrix(convection->get_concentration_matrix(), el_distribution, el_4_loc);

	        //Supresses possibility to combine reactions
	        /*Semchem_reactions = NULL;
	        sorptions = NULL;*/
		} else
	    if (reactions_it->type() == Pade_approximant::input_type) {
            decayRad = new Pade_approximant(init_mesh, *reactions_it, subst_names_ );
	        convection->get_par_info(el_4_loc, el_distribution);
	        decayRad->set_dual_porosity(convection->get_dual_porosity());
	        static_cast<Pade_approximant *> (decayRad) -> modify_reaction_matrix();
	        decayRad->set_concentration_matrix(convection->get_concentration_matrix(), el_distribution, el_4_loc);

	        //Supresses possibility to combine reactions
	        /*Semchem_reactions = NULL;
	        sorptions = NULL;*/
	    } else
	    if (reactions_it->type() == Semchem_interface::input_type ) {
	        Semchem_reactions = new Semchem_interface(0.0, mesh_, n_subst_, convection->get_dual_porosity()); //(mesh->n_elements(),convection->get_concentration_matrix(), mesh);
	        Semchem_reactions->set_el_4_loc(el_4_loc);
	        Semchem_reactions->set_concentration_matrix(convection->get_concentration_matrix(), el_distribution, el_4_loc);

	        /*decayRad = NULL;
	        sorptions = NULL;*/
	    } else {
	        xprintf(UsrErr, "Wrong reaction type.\n");
	    }
	} else {
	    decayRad = NULL;
	    Semchem_reactions = NULL;
	}

	Input::Iterator<Input::Record> sorptions_it = in_rec.find<Input::Record>("adsorptions");
	if (sorptions_it){
        convection->get_par_info(el_4_loc, el_distribution);
        // Part for mobile zone description follows.
	    sorptions = new Sorption(init_mesh, *sorptions_it, subst_names_);
	    sorptions->set_dual_porosity(convection->get_dual_porosity());
	    sorptions->set_porosity(&(convection->get_data()->por_m), &(convection->get_data()->por_imm)); //, &(convection->get_data()->por_imm));
	    sorptions->set_phi(&(convection->get_data()->phi));
	    sorptions->prepare_inputs(*sorptions_it);
	    double ***conc_matrix = convection->get_concentration_matrix();
	    sorptions->set_concentration_matrix(conc_matrix[MOBILE], el_distribution, el_4_loc);
	    sorptions->set_sorb_conc_array(el_distribution->lsize());

	    if(convection->get_dual_porosity()){
	    	sorptions_immob = new Sorption_dp(init_mesh, *sorptions_it, subst_names_);
		    sorptions_immob->set_dual_porosity(convection->get_dual_porosity());
	    	sorptions_immob->set_nr_transp(n_subst_);
	    	sorptions_immob->set_porosity(&(convection->get_data()->por_m), &(convection->get_data()->por_imm));
	    	sorptions_immob->set_phi(&(convection->get_data()->phi));
		    sorptions_immob->prepare_inputs(*sorptions_it);
		    sorptions_immob->set_concentration_matrix(conc_matrix[MOBILE], el_distribution, el_4_loc);
		    sorptions_immob->set_immob_concentration_matrix(conc_matrix[IMMOBILE], el_distribution, el_4_loc);
		    sorptions_immob->set_sorb_conc_array(el_distribution->lsize());
	    }
	  } else{
	    sorptions = NULL;
	    sorptions_immob = NULL;
	}
	
	output_mark_type = convection->mark_type() | TimeGovernor::marks().type_fixed_time() | TimeGovernor::marks().type_output();
    time_ = new TimeGovernor(in_rec.val<Input::Record>("time"), output_mark_type );

}

TransportOperatorSplitting::~TransportOperatorSplitting()
{
    //delete field_output;
    delete convection;
    if (decayRad) delete decayRad;
    if (sorptions) delete sorptions;
    if (sorptions_immob) delete sorptions_immob;
    if (Semchem_reactions) delete Semchem_reactions;
    delete time_;
}




void TransportOperatorSplitting::output_data(){

    if (time_->is_current(output_mark_type)) {
        
        START_TIMER("TOS-output data");
        DBGMSG("\nTOS: output time: %f\n", time_->t());
        
        convection->output_data();
    }
}



void TransportOperatorSplitting::update_solution() {

	static bool first_time_call = true;

	if (first_time_call)
	{
		convection->mass_balance();
		first_time_call = false;
	}

    time_->next_time();
    time_->view("TOS");    //show time governor
    
    convection->set_target_time(time_->t());
	if (decayRad) decayRad->set_time_step(convection->time().estimate_dt());
	if (sorptions) sorptions->set_time_step(convection->time().estimate_dt());
	if (sorptions_immob && convection->get_dual_porosity()) sorptions_immob->set_time_step(convection->time().estimate_dt());
	// TODO: update Semchem time step here!!
	if (Semchem_reactions) Semchem_reactions->set_timestep(convection->time().estimate_dt());

        
    xprintf( Msg, "TOS: time: %f        CONVECTION: time: %f      dt_estimate: %f\n", 
             time_->t(), convection->time().t(), convection->time().estimate_dt() );
    
    START_TIMER("TOS-one step");
    int steps=0;
    while ( convection->time().lt(time_->t()) )
    {
        steps++;
	    // one internal step
	    convection->compute_one_step();
		if (sorptions_immob && convection->get_dual_porosity()) sorptions_immob->transport_dual_porosity(); // belongs to completely different place
	    if(decayRad) decayRad->compute_one_step();
	    if(Semchem_reactions) Semchem_reactions->compute_one_step();
	    if(sorptions) sorptions->compute_one_step();//equilibrial sorption at the end of simulated time-step
	    if(sorptions_immob && convection->get_dual_porosity()) sorptions_immob->compute_one_step();
	}
    END_TIMER("TOS-one step");
    
    xprintf( Msg, "CONVECTION: steps: %d\n",steps);
}





void TransportOperatorSplitting::set_velocity_field(const MH_DofHandler &dh)
{
    mh_dh = &dh;
	convection->set_velocity_field( dh );
};



void TransportOperatorSplitting::get_parallel_solution_vector(Vec &vec){
	convection->compute_one_step();
};



void TransportOperatorSplitting::get_solution_vector(double * &x, unsigned int &a){
	convection->compute_one_step();
};



void TransportOperatorSplitting::set_eq_data(Field< 3, FieldValue<3>::Scalar >* cross_section)
{
    convection->set_cross_section_field(cross_section);

    /*if (Semchem_reactions != NULL) {
        Semchem_reactions->set_cross_section(cross_section);
        Semchem_reactions->set_sorption_fields(&convection->get_data()->por_m, &convection->get_data()->por_imm, &convection->get_data()->phi);
    }
  if (sorptions != NULL)
  {
	  sorptions->set_porosity(&(convection->get_data()->por_m),&(convection->get_data()->por_imm));
  }*/
}



void TransportOperatorSplitting::calc_fluxes(vector<vector<double> > &bcd_balance, vector<vector<double> > &bcd_plus_balance, vector<vector<double> > &bcd_minus_balance)
{
    convection->calc_fluxes(bcd_balance, bcd_plus_balance, bcd_minus_balance);
}

void TransportOperatorSplitting::calc_elem_sources(vector<vector<double> > &mass, vector<vector<double> > &src_balance)
{
    convection->calc_elem_sources(mass, src_balance);
}





