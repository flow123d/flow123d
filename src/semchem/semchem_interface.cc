//---------------------------------------------------------------------------

#include "reaction/reaction.hh"

#include "system/system.hh"
#include "semchem/che_semchem.h"
#include "semchem/semchem_interface.hh"
#include "transport/transport.h"
#include "mesh/mesh.h"

using namespace std;

//---------------------------------------------------------------------------
//  GLOBALNI PROMENNE
//---------------------------------------------------------------------------
struct TS_prm	G_prm;
struct TS_lat 	*P_lat;
struct TS_che	*P_che;

//---------------------------------------------------------------------------

Input::Type::Record & Specie::get_input_type()
{
	using namespace Input::Type;
	static Record rec("Isotope", "Definition of information about a single isotope.");

	if (!rec.is_finished()) {
		rec.declare_key("identifier", Integer(), Default::obligatory(),
						"Identifier of the isotope.");
		rec.declare_key("half_life", Double(), Default::obligatory(),
						"Half life parameter.");
		/*rec.declare_key("next", Array(Integer()), Default(0),
						"Identifiers of childern in decay chain.");
		rec.declare_key("bifurcation", Array(Double), Default(0),
						"Fractions of division decay chain into branches.");
		rec.declare_key("kinetic constant", Double(), Default(1.0),
						"Kinetic conxtant appropriate to described first order reaction.");*/

		rec.finish();
	}
	return rec;
}

Input::Type::Record & General_reaction::get_input_type()
{
	using namespace Input::Type;
	static Record rec("Isotope", "Definition of information about a single isotope.");

	if (!rec.is_finished()) {
	    rec.derive_from(Reaction::get_input_type());

        //rec.declare_key("general_reaction", Array( Linear_reaction::get_one_decay_substep() ), Default::optional(),
        //        "Description of general chemical reactions.");

		rec.declare_key("identifier", Integer(), Default::obligatory(),
						"Identifier of the isotope.");
		rec.declare_key("half_life", Double(), Default::obligatory(),
						"Half life parameter.");
		/*rec.declare_key("next", Array(Integer()), Default(0),
						"Identifiers of childern in decay chain.");
		rec.declare_key("bifurcation", Array(Double), Default(0),
						"Fractions of division decay chain into branches.");
		rec.declare_key("kinetic constant", Double(), Default(1.0),
						"Kinetic conxtant appropriate to described first order reaction.");*/

		rec.finish();
	}
	return rec;
}

Input::Type::AbstractRecord & Semchem_interface::get_input_type()
{
	using namespace Input::Type;
	static AbstractRecord rec("Semchem_module", "Declares infos valid for all reactions.");

	if (!rec.is_finished()) {
		rec.declare_key("precision", Integer(), Default::obligatory(), //(1),
						"How accurate should the simulation be, decimal places(?).");
		rec.declare_key("temperature", Double(), Default::obligatory(), //(298.0),
						"Isothermal reaction, thermodynamic temperature.");
		rec.declare_key("temp_Gf", Double(), Default::obligatory(), //(298.0),
						"Thermodynamic parameter.");
		rec.declare_key("param_Afi", Double(), Default::obligatory(), //(0.391475),
						"Thermodynamic parameter.");
		rec.declare_key("param_b", Double(), Default::obligatory(), //(1.2),
						"Thermodynamic parameter.");
		rec.declare_key("epsilon", Double(), Default::obligatory(), //(1.2),
						"Thermodynamic parameter.");
		rec.declare_key("time_steps", Integer(), Default::obligatory(), //(10),
						"Simulation parameter.");
		rec.declare_key("slow_kinetic_steps", Integer(), Default::obligatory(), //(1),
						"Simulation parameter.");

		rec.finish();

		//TransportOperatorSplitting::get_input_type();
		//TransportDG::get_input_type();

		rec.no_more_descendants();
	}
	return rec;
}

Semchem_interface::Semchem_interface(double timeStep, Mesh * mesh, int nrOfSpecies, bool dualPorosity)
	:semchem_on(false), dual_porosity_on(false), mesh_(NULL), fw_chem(NULL)
{

  //temporary semchem output file name
  std::string semchem_output_fname = IONameHandler::get_instance()->get_output_file_name("./semchem_output.out");
  xprintf(Msg,"Semchem output file name is %s\n",semchem_output_fname.c_str());

  //char *Semchem_output_file;
  //= semchem_output_fname.c_str();
  //Semchem_output_file = (char *)xmalloc(sizeof(char)*(semchem_output_fname.length() + 1));
  //sprintf(Semchem_output_file,"%s",semchem_output_fname.c_str());
  //strcpy(Semchem_output_file,semchem_output_fname.c_str());
  //Semchem_output_file[semchem_output_fname.length()] = "\0";

  this->set_fw_chem(semchem_output_fname); //DOES NOT WORK ((const char*)semchem_output_fname.c_str());
  this->set_chemistry_computation();
  if(semchem_on == true) ctiich();
  set_dual_porosity();
  set_mesh_(mesh);
  set_nr_of_elements(mesh_->n_elements());
  return;
}

/*Semchem_interface::~Semchem_interface(void)
{
	if(fw_chem != NULL)
	{
		free(fw_chem);
		fw_chem = NULL;
	}
	return;
}*/

void Semchem_interface::compute_one_step(void)
{
	if(semchem_on == true)
	{
		for (int loc_el = 0; loc_el < distribution->lsize(); loc_el++)
		{
			START_TIMER("semchem_step");
	   	   this->compute_reaction(dual_porosity_on, mesh_->element(el_4_loc[loc_el]), loc_el, concentration_matrix);
	   	   END_TIMER("semchem_step");
		}
	}
}

//---------------------------------------------------------------------------
//                 FUNKCE NA VYPOCET CHEMIE PRO ELEMENT V TRANSPORTU
//---------------------------------------------------------------------------
//void Semchem_interface::compute_reaction(bool porTyp, ElementIter ppelm, int poradi, double **conc_mob_arr, double **conc_immob_arr)
void Semchem_interface::compute_reaction(bool porTyp, ElementIter ppelm, int poradi, double ***conc)
{
  FILE *fw, *stream, *fr;
   int i, j; //, poradi;
   int krok;
   int poc_krok;
   double celkova_molalita;
   char *vystupni_soubor; //[] = "./output/semchem_output.txt";
   double **conc_mob_arr = conc[MOBILE];
   double **conc_immob_arr = conc[IMMOBILE];
   double pomoc, n;

   vystupni_soubor = fw_chem;
   //==================================================================
   // ----------- ALOKACE POLE PRO KONCENTRACE Z FLOWA ----------------
   //==================================================================

   ASSERT(P_lat != NULL,"\nP_lat NENI ALOKOVANE\n");

  //==================================================================
   // ----------------- NEJPRVE PRO MOBILNI PORY ----------------------
   //==================================================================
   n = (ppelm->material->por_m) / (1 - ppelm->material->por_m); //asi S/V jako zze splocha

   switch (ppelm->dim) { //objem se snad na nic nepouzzi:va:
	  case 1 :
	  case 2 :
	  case 3 : pomoc = (ppelm->volume) * (ppelm->material->por_m); break;
	default:
	  pomoc = 1.0;
   }

   G_prm.objem = pomoc; //objem * mobilni: porozita
   G_prm.splocha = (pomoc / ppelm->material->por_m) * (ppelm->material->phi) * (1 - ppelm->material->por_m - ppelm->material->por_imm);
   celkova_molalita=0.0;
   poc_krok=1;

    //------------PREDANI VSECH INFORMACI Z FLOWA DO SEMCHEMU----------
   for ( i=0 ; i<(G_prm.pocet_latekvefazi); i++)
   {
	 P_lat[i].m0 = (double)((conc_mob_arr[i][poradi])) / (P_lat[i].M);
	 celkova_molalita += (P_lat[i].m0);
   }
   G_prm.deltaT = time_step/G_prm.cas_kroku; // dosazeni "spravneho" casoveho kroku

    //-----------------------VYPOCET CHEMIE----------------------------
   if (celkova_molalita > 1e-16)
   {
	  for (krok = 1; krok <= G_prm.cas_kroku; krok++)
	  {
		 che_pocitej_soubor(vystupni_soubor, &poc_krok);
		 for (i = 0; i < G_prm.pocet_latekvefazi; i++) {
		   P_lat[i].m0 = P_lat[i].m;
		 }
		 che_presun_poc_p_();
	  }

    //------------PREDANI VSECH INFORMACI ZE SEMCHEMU DO FLOWA---------
	  for ( i=0 ; i<G_prm.pocet_latekvefazi; i++)
	  {
			conc_mob_arr[i][poradi] = (double)(P_lat[i].m0 * P_lat[i].M);
	  }
	}

    //==================================================================
    // ----------------- POTE PRO IMOBILNI PORY ------------------------
    //==================================================================
   if(porTyp == true) {
   switch (ppelm->dim) { //objem se snad na nic nepouzzi:va:
	  case 1 :
	  case 2 :
	  case 3 : pomoc = ppelm->volume * (ppelm->material->por_imm); break;
	default:
	  pomoc = 1.0;
    }
   G_prm.objem = pomoc;
   G_prm.splocha = (pomoc / ppelm->material->por_imm) * (1 - ppelm->material->phi) * (1 - ppelm->material->por_m - ppelm->material->por_imm);
   celkova_molalita = 0.0;
   poc_krok=1;
   
    //------------PREDANI VSECH INFORMACI Z FLOWA DO SEMCHEMU----------
   for ( i=0 ; i<G_prm.pocet_latekvefazi ; i++)
   {
	 P_lat[i].m0 = (double)(conc_immob_arr[i][poradi] / P_lat[i].M);
	 celkova_molalita += P_lat[i].m0;
   }

    //-----------------------VYPOCET CHEMIE----------------------------
  if (celkova_molalita > 1e-16)
   {
	  //if (G_prm.vypisy==1) che_nadpis__soubor(fw_chem);

	  for (krok = 1; krok<=G_prm.cas_kroku; krok++)
	  {
		 che_pocitej_soubor(vystupni_soubor,&poc_krok);
		 che_presun_poc_p_();
	  }

    //------------PREDANI VSECH INFORMACI ZE SEMCHEMU DO FLOWA---------
	  for ( i=0 ; i<G_prm.pocet_latekvefazi ; i++)
	  {
			conc_immob_arr[i][poradi] = (P_lat[i].m0 * P_lat[i].M); //m nebo m0, co?
	  }
     }
    }else{
	  ;
    }
}

void Semchem_interface::set_timestep(double new_timestep)
{
	this->time_step = new_timestep;
	return;
}

void Semchem_interface::set_chemistry_computation(void)
{
	this->semchem_on = OptGetBool("Semchem_module", "Compute_reactions", "no");
	return;
}

void Semchem_interface::set_dual_porosity()
{
	this->dual_porosity_on = OptGetBool("Transport", "Dual_porosity", "no");
	return;
}

void Semchem_interface::set_nr_of_elements(int nrOfElements)
{
	this->nr_of_elements = nrOfElements;
	return;
}

void Semchem_interface::set_concentration_matrix(double ***ConcentrationMatrix, Distribution *conc_distr, int *el_4_loc)
{
	concentration_matrix = ConcentrationMatrix;
	distribution = conc_distr;
	return;
}

void Semchem_interface::set_el_4_loc(int *el_for_loc)
{
	el_4_loc = el_for_loc;
	return;
}

void Semchem_interface::set_mesh_(Mesh *mesh)
{
	mesh_ = mesh;
	return;
}

void Semchem_interface::set_fw_chem(std::string semchem_output_file) //(const char* semchem_output_file)
{
	//fw_chem = (char*)xmalloc(sizeof((semchem_output_file)/sizeof(char)+1)*sizeof(char));
	fw_chem = (char*)xmalloc(semchem_output_file.length()+1);
	strcpy(fw_chem,semchem_output_file.c_str());
	xprintf(Msg,"Output file for Semchem is %s\n",fw_chem);
	return;
}
