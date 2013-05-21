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
#include "reaction/isotherm.hh"
#include "reaction/sorption.hh"

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
//	.declare_key("output_stream", OutputTime::input_type, Default::obligatory(),
//			"Parameters of output stream.")
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
    .declare_key("bc_data", Array(TransportOperatorSplitting::EqData().boundary_input_type()
    		.declare_key("old_boundary_file", IT::FileName::input(), "Input file with boundary conditions (obsolete).")
    		.declare_key("bc_times", Array(Double()), Default::optional(),
    				"Times for changing the boundary conditions (obsolete).")
    		), IT::Default::obligatory(), "")
    .declare_key("bulk_data", Array(TransportOperatorSplitting::EqData().bulk_input_type()),
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
    const int nsubst = n_substances();
    vector<vector<double> > bcd_balance(nsubst),
    		bcd_plus_balance(nsubst),
    		bcd_minus_balance(nsubst),
    		mass(nsubst),
    		src_balance(nsubst);


    for (int i=0; i<nsubst; i++)
    {
    	bcd_balance[i].resize(mesh_->region_db().boundary_size());
    	bcd_plus_balance[i].resize(mesh_->region_db().boundary_size());
    	bcd_minus_balance[i].resize(mesh_->region_db().boundary_size());
    	for (unsigned int j=0; j<mesh_->region_db().boundary_size(); j++)
    	{
    		bcd_balance[i][j] = 0;
    		bcd_plus_balance[i][j] = 0;
    		bcd_minus_balance[i][j] = 0;
    	}

    	mass[i].resize(mesh_->region_db().bulk_size());
    	src_balance[i].resize( mesh_->region_db().bulk_size());
		for (unsigned int j=0; j<src_balance[i].size(); j++)
		{
			mass[i][j] = 0;
			src_balance[i][j] = 0;
		}
	}
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
    string bc_head_format = "# %-*s%-*s%-*s%-*s%-*s%-*s%-*s\n",
           bc_format = "%*s%-*d%-*s%-*s  %-*g%-*g%-*g%-*g\n",
           bc_total_format = "# %-*s%-*s%-*g%-*g%-*g\n";
    s << setw((w*c+wl-14)/2) << setfill('-') << "--"; //drawing half line
    fprintf(balance_output_file,"# %s MASS BALANCE %s\n",s.str().c_str(), s.str().c_str());
    fprintf(balance_output_file,"# Time of computed mass balance: %f\n\n\n",time().t());

    //BOUNDARY
    fprintf(balance_output_file,"# Boundary mass balance:\n");
    fprintf(balance_output_file,bc_head_format.c_str(),w,"[boundary_id]",wl,"[label]",
                            w,"[substance]",w,"[total_balance]",w,"[total_outflow]",w,"[total_inflow]",w,"[time]");
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
					w, bcd_minus_balance[sbi][reg->boundary_idx()], w, time().t());
    	}
    }
    //total boundary balance
    fprintf(balance_output_file,"# %s\n",s.str().c_str());  // drawing long line
    for (int sbi=0; sbi<nsubst; sbi++)
    	fprintf(balance_output_file, bc_total_format.c_str(),w+wl,"total boundary balance",
                w+2,substance_names()[sbi].c_str(),w,bcd_total_balance[sbi], w, bcd_total_outflow[sbi], w, bcd_total_inflow[sbi]);
    fprintf(balance_output_file, "\n\n");


    //SOURCES
    string src_head_format = "# %-*s%-*s%-*s%-*s%-*s%-*s%-*s\n",
           src_format = "%*s%-*d%-*s%-*s  %-*g%-*g%-*s%-*g\n",
           src_total_format = "# %-*s%-*s%-*g%-*g\n";
    //computing water balance of sources
    fprintf(balance_output_file,"# Source fluxes over material subdomains:\n");   //head
    fprintf(balance_output_file,src_head_format.c_str(),w,"[region_id]",wl,"[label]",
                            w,"[substance]",w,"[total_mass]",w,"[total_sources]",w,"",w,"[time]");
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
					w,src_balance[sbi][reg->bulk_idx()],
					w,"", w,time().t());
    	}
    }
    //total sources balance
    fprintf(balance_output_file,"# %s\n",s.str().c_str());  //drawing long line
    for (int sbi=0; sbi<nsubst; sbi++)
    	fprintf(balance_output_file, src_total_format.c_str(),w+wl,"total sources balance",
                w+2,substance_names()[sbi].c_str(),
                w,mass_total[sbi],
                w,src_total_balance[sbi]);
    fprintf(balance_output_file, "\n\n");
}




RegionSet TransportOperatorSplitting::EqData::read_boundary_list_item(Input::Record rec) {
	// Base method EqDataBase::read_boundary_list_item must be called first!
	RegionSet domain = EqDataBase::read_boundary_list_item(rec);
    FilePath bcd_file;

    // read transport boundary conditions using old file format .tbc
    if (rec.opt_val("old_boundary_file", bcd_file) )
        OldBcdInput::instance()->read_transport(bcd_file, bc_conc);

    return domain;
}




TransportOperatorSplitting::EqData::EqData() : TransportEqData("TransportOperatorSplitting")
{
	ADD_FIELD(por_imm, "Immobile porosity", Default("0"));
	ADD_FIELD(alpha, "Coefficients of non-equilibrium exchange.", Default("0"));
	ADD_FIELD(sorp_type, "Type of sorption.", Default("1"));
	ADD_FIELD(sorp_coef0, "Coefficient of sorption.", Default("0"));
	ADD_FIELD(sorp_coef1, "Coefficient of sorption.", Default("0"));
	ADD_FIELD(phi, "Solid / solid mobile.", Default("0.5"));

}


TransportOperatorSplitting::TransportOperatorSplitting(Mesh &init_mesh, const Input::Record &in_rec)
: TransportBase(init_mesh, in_rec)
{
	Distribution *el_distribution;
	int *el_4_loc;

    // double problem_save_step = OptGetDbl("Global", "Save_step", "1.0");

	convection = new ConvectionTransport(*mesh_, data, in_rec);

	Input::Iterator<Input::AbstractRecord> reactions_it = in_rec.find<Input::AbstractRecord>("reactions");
	if ( reactions_it ) {
		if (reactions_it->type() == Linear_reaction::input_type ) {
	        decayRad =  new Linear_reaction(init_mesh, *reactions_it, convection->get_substance_names());
	        convection->get_par_info(el_4_loc, el_distribution);
	        decayRad->set_dual_porosity(convection->get_dual_porosity());
	        static_cast<Linear_reaction *> (decayRad) -> modify_reaction_matrix();
	        decayRad->set_concentration_matrix(convection->get_prev_concentration_matrix(), el_distribution, el_4_loc);

	        Semchem_reactions = NULL;
	        sorptions = NULL;
		} else
	    if (reactions_it->type() == Pade_approximant::input_type ) {
                decayRad = new Pade_approximant(init_mesh, *reactions_it, convection->get_substance_names());
	        convection->get_par_info(el_4_loc, el_distribution);
	        decayRad->set_dual_porosity(convection->get_dual_porosity());
	        static_cast<Pade_approximant *> (decayRad) -> modify_reaction_matrix();
	        decayRad->set_concentration_matrix(convection->get_prev_concentration_matrix(), el_distribution, el_4_loc);

	        Semchem_reactions = NULL;
	        sorptions = NULL;
	    } else
	    if (reactions_it->type() == Semchem_interface::input_type ) {
	        Semchem_reactions = new Semchem_interface(0.0, mesh_, convection->get_n_substances(), convection->get_dual_porosity()); //(mesh->n_elements(),convection->get_concentration_matrix(), mesh);
	        Semchem_reactions->set_el_4_loc(el_4_loc);
	        Semchem_reactions->set_concentration_matrix(convection->get_prev_concentration_matrix(), el_distribution, el_4_loc);

	        decayRad = NULL;
	        sorptions = NULL;
	    } else
	    if (reactions_it->type() == Sorption::input_type ){
	    	sorptions = new Sorption(init_mesh, *reactions_it, convection->get_substance_names());
	        convection->get_par_info(el_4_loc, el_distribution);
	        sorptions->set_dual_porosity(convection->get_dual_porosity());

	        sorptions->set_concentration_matrix(convection->get_prev_concentration_matrix(), el_distribution, el_4_loc);

	        double** sorb_conc_array;
	        sorb_conc_array = (double**) xmalloc(convection->get_n_substances() * sizeof(double*));
	        for (int sbi = 0; sbi < convection->get_n_substances(); sbi++)
	        {
	          sorb_conc_array[sbi] = (double*) xmalloc(el_distribution->lsize() * sizeof(double));
	          for (int i = 0; i < el_distribution->lsize(); i++)
	          {
	            sorb_conc_array[sbi][i] = 0.0;
	          }
	        }
	        sorptions->set_sorb_conc_array(sorb_conc_array);

	        decayRad = NULL;
	        Semchem_reactions = NULL;
	    } else{
	        xprintf(UsrErr, "Wrong reaction type.\n");
	    }
	} else {
	    decayRad = NULL;
	    Semchem_reactions = NULL;
	    sorptions = NULL;
	}
	
        time_ = new TimeGovernor(in_rec.val<Input::Record>("time"), this->mark_type());
        output_mark_type = this->mark_type() | time_->marks().type_fixed_time() | time_->marks().type_output();

        time_->marks().add_time_marks(0.0,
            in_rec.val<Input::Record>("output").val<double>("save_step"),
            time_->end_time(), output_mark_type );
	// TODO: this has to be set after construction of transport matrix !!


	// register output vectors from convection
	double ***out_conc = convection->get_out_conc();
	vector<string> substance_name = convection->get_substance_names();

	// TODO: Add corresponding record to the in_rec
	Input::Record output_rec = in_rec.val<Input::Record>("output");

	//field_output = new OutputTime(mesh_, output_rec.val<Input::Record>("output_stream"));
	//field_output = OutputTime::output_stream(output_rec.val<Input::Record>("output_stream"));

    for(int subst_id=0; subst_id < convection->get_n_substances(); subst_id++) {
         // TODO: What about output also other "phases", IMMOBILE and so on.
         std::string subst_name = substance_name[subst_id] + "_mobile";
         double *data = out_conc[MOBILE][subst_id];
         this->output_streams[data] = OutputTime::register_elem_data<double>(mesh_, subst_name, "M/L^3",
                 output_rec, data , mesh_->n_elements(), time_->t());
    }
    // write initial condition
    convection->output_vector_gather();
    //if(field_output) field_output->write_data(time_->t());

}

TransportOperatorSplitting::~TransportOperatorSplitting()
{
    //delete field_output;
    delete convection;
    if (decayRad) delete decayRad;
    if (sorptions) delete sorptions;
    if (Semchem_reactions) delete Semchem_reactions;
    delete time_;
}




void TransportOperatorSplitting::output_data(){

    if (time_->is_current(output_mark_type)) {
        
        START_TIMER("TOS-output data");
        
        DBGMSG("\nTOS: output time: %f\n", time_->t());
        convection->output_vector_gather();
        //if(field_output) field_output->write_data(time_->t());
        for(std::map<void*, OutputTime*>::iterator it = this->output_streams.begin();
                it != this->output_streams.end();
                ++it)
        {
            ((OutputTime*)it->second)->set_data_time(it->first, time_->t());
        }
        mass_balance();
        
        //for synchronization when measuring time by Profiler
        MPI_Barrier(MPI_COMM_WORLD);
    }
}



void TransportOperatorSplitting::update_solution() {


    time_->next_time();
    //time_->view("TOS");    //show time governor
    
    convection->set_target_time(time_->t());

	if (decayRad) decayRad->set_time_step(convection->time().estimate_dt());
	if (sorptions) sorptions->set_time_step(convection->time().estimate_dt());
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
	    // Calling linear reactions and Semchem, temporarily commented
	    if(decayRad) decayRad->compute_one_step();
	    if(Semchem_reactions) Semchem_reactions->compute_one_step();
	    if(sorptions) sorptions->compute_one_step();//equilibrial sorption at the end of simulated time-step
	}
    END_TIMER("TOS-one step");
    
    xprintf( Msg, "CONVECTION: steps: %d\n",steps);
}






void TransportOperatorSplitting::set_velocity_field(const MH_DofHandler &dh)
{
    mh_dh = &dh;
	convection->set_flow_field_vector( dh );
};


void TransportOperatorSplitting::get_parallel_solution_vector(Vec &vec){
	convection->compute_one_step();
};

void TransportOperatorSplitting::get_solution_vector(double * &x, unsigned int &a){
	convection->compute_one_step();
};

void TransportOperatorSplitting::set_eq_data(Field< 3, FieldValue<3>::Scalar >* cross_section)
{
  data.cross_section = cross_section;
  if (Semchem_reactions != NULL) {
	  Semchem_reactions->set_cross_section(cross_section);
	  Semchem_reactions->set_sorption_fields(&data.por_m, &data.por_imm, &data.phi);
  }
}

unsigned int TransportOperatorSplitting::n_substances()
{
	return convection->get_n_substances();
}

vector<string> &TransportOperatorSplitting::substance_names()
{
	return convection->get_substance_names();
}


void TransportOperatorSplitting::calc_fluxes(vector<vector<double> > &bcd_balance, vector<vector<double> > &bcd_plus_balance, vector<vector<double> > &bcd_minus_balance)
{
	double ***solution = convection->get_conc();
	int *el_4_loc, *row_4_el;
	Distribution *el_ds;
	double mass_flux[n_substances()];

	convection->get_par_info(el_4_loc, el_ds);
	row_4_el = convection->get_row_4_el();

    FOR_BOUNDARIES(mesh_, bcd) {

        // !! there can be more sides per one boundary
    	int index = row_4_el[bcd->side()->element().index()];
    	if (!el_ds->is_local(index)) continue;

		double water_flux = mh_dh->side_flux(*(bcd->side()));
		for (unsigned int sbi=0; sbi<n_substances(); sbi++)
			mass_flux[sbi] = water_flux*solution[MOBILE][sbi][index-el_ds->begin()];

        Region r = bcd->region();
        if (! r.is_valid()) xprintf(Msg, "Invalid region, ele % d, edg: % d\n", bcd->bc_ele_idx_, bcd->edge_idx_);
        unsigned int bc_region_idx = r.boundary_idx();

        for (unsigned int sbi=0; sbi<n_substances(); sbi++)
        {
        	bcd_balance[sbi][bc_region_idx] += mass_flux[sbi];

        	if (mass_flux[sbi] > 0) bcd_plus_balance[sbi][bc_region_idx] += mass_flux[sbi];
        	else bcd_minus_balance[sbi][bc_region_idx] += mass_flux[sbi];
        }
    }

}

void TransportOperatorSplitting::calc_elem_sources(vector<vector<double> > &mass, vector<vector<double> > &src_balance)
{
	int *el_4_loc, *row_4_el;
	Distribution *el_ds;
	double ***solution = convection->get_conc();

	convection->get_par_info(el_4_loc, el_ds);
	row_4_el = convection->get_row_4_el();

	for (unsigned int sbi=0; sbi<n_substances(); sbi++)
	{
		double *sources = convection->get_sources(sbi);

		FOR_ELEMENTS(mesh_,elem)
		{
			int index = row_4_el[elem.index()];
			if (el_ds->is_local(index))
			{
				mass[sbi][elem->element_accessor().region().bulk_idx()] += solution[MOBILE][sbi][index-el_ds->begin()]*elem->measure();
				src_balance[sbi][elem->element_accessor().region().bulk_idx()] += sources[index-el_ds->begin()]*elem->measure();
			}
		}
	}
}

