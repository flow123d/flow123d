//---------------------------------------------------------------------------

#include "reaction/reaction.hh"

#include "system/system.hh"
#include "system/sys_profiler.hh"
#include "io/read_ini.h"

#include "semchem/che_semchem.h"
#include "semchem/semchem_interface.hh"
//#include "transport/transport.h"
#include "mesh/mesh.h"
#include "fields/field_base.hh"
#include "fields/field_values.hh"

#define MOBILE 0
#define IMMOBILE 1

using namespace std;

//---------------------------------------------------------------------------
//  GLOBALNI PROMENNE
//---------------------------------------------------------------------------
struct TS_prm	G_prm;
struct TS_lat 	*P_lat;
struct TS_che	*P_che;

//---------------------------------------------------------------------------
namespace it = Input::Type;

it::Record Specie::input_type = it::Record("Isotope", "Definition of information about a single isotope.")
	.declare_key("identifier", it::Integer(), it::Default::obligatory(),
						"Identifier of the isotope.")
	.declare_key("half_life", it::Double(), it::Default::obligatory(),
						"Half life parameter.");
		/*rec.declare_key("next", Array(Integer()), Default(0),
						"Identifiers of childern in decay chain.");
		rec.declare_key("bifurcation", Array(Double), Default(0),
						"Fractions of division decay chain into branches.");
		rec.declare_key("kinetic constant", Double(), Default(1.0),
						"Kinetic conxtant appropriate to described first order reaction.");*/


it::Record General_reaction::input_type = it::Record("Isotope", "Definition of information about a single isotope.")
	.derive_from(ReactionTerm::input_type)
        //rec.declare_key("general_reaction", Array( Linear_reaction::get_one_decay_substep() ), Default::optional(),
        //        "Description of general chemical reactions.");
	.declare_key("identifier", it::Integer(), it::Default::obligatory(),
						"Identifier of the isotope.")
	.declare_key("half_life", it::Double(), it::Default::obligatory(),
						"Half life parameter.");
		/*rec.declare_key("next", Array(Integer()), Default(0),
						"Identifiers of childern in decay chain.");
		rec.declare_key("bifurcation", Array(Double), Default(0),
						"Fractions of division decay chain into branches.");
		rec.declare_key("kinetic constant", Double(), Default(1.0),
						"Kinetic conxtant appropriate to described first order reaction.");*/


it::AbstractRecord Semchem_interface::input_type = it::AbstractRecord("Semchem_module", "Declares infos valid for all reactions.")
	.declare_key("precision", it::Integer(), it::Default::obligatory(), //(1),
						"How accurate should the simulation be, decimal places(?).")
	.declare_key("temperature", it::Double(), it::Default::obligatory(), //(298.0),
						"Isothermal reaction, thermodynamic temperature.")
	.declare_key("temp_gf", it::Double(), it::Default::obligatory(), //(298.0),
						"Thermodynamic parameter.")
	.declare_key("param_afi", it::Double(), it::Default::obligatory(), //(0.391475),
						"Thermodynamic parameter.")
	.declare_key("param_b", it::Double(), it::Default::obligatory(), //(1.2),
						"Thermodynamic parameter.")
	.declare_key("epsilon", it::Double(), it::Default::obligatory(), //(1.2),
						"Thermodynamic parameter.")
	.declare_key("time_steps", it::Integer(), it::Default::obligatory(), //(10),
						"Simulation parameter.")
	.declare_key("slow_kinetic_steps", it::Integer(), it::Default::obligatory(), //(1),
						"Simulation parameter.");


Semchem_interface::Semchem_interface(double timeStep, Mesh * mesh, int nrOfSpecies, bool dualPorosity)
	:semchem_on(false), dual_porosity_on(false), fw_chem(NULL), mesh_(NULL), cross_section(cross_section)
{

  //temporary semchem output file name
  std::string semchem_output_fname = FilePath("semchem_output.out", FilePath::output_file);
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

void Semchem_interface::set_cross_section(Field< 3 , FieldValue< 3  >::Scalar >* cross_section)
{
  this->cross_section = cross_section;
}

void Semchem_interface::set_sorption_fields(Field<3, FieldValue<3>::Scalar> *por_m_,
		Field<3, FieldValue<3>::Scalar> *por_imm_,
		Field<3, FieldValue<3>::Scalar> *phi_)
{
	por_m = por_m_;
	por_imm = por_imm_;
	phi = phi_;
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

void Semchem_interface::update_solution(void)
{
	if(semchem_on == true)
	{
		for (unsigned int loc_el = 0; loc_el < distribution->lsize(); loc_el++)
		{
			START_TIMER("semchem-one step");
	   	   this->compute_reaction(dual_porosity_on, mesh_->element(el_4_loc[loc_el]), loc_el, concentration_matrix);
	   	   END_TIMER("semchem-one step");
		}
	}
}

//---------------------------------------------------------------------------
//                 FUNKCE NA VYPOCET CHEMIE PRO ELEMENT V TRANSPORTU
//---------------------------------------------------------------------------
//void Semchem_interface::compute_reaction(bool porTyp, ElementIter ppelm, int poradi, double **conc_mob_arr, double **conc_immob_arr)
void Semchem_interface::compute_reaction(bool porTyp, ElementIter ppelm, int poradi, double ***conc)
{
   //FILE *fw, *stream, *fr;
   int i; //, j, poradi;
   int krok;
   int poc_krok;
   double celkova_molalita;
   char *vystupni_soubor; //[] = "./output/semchem_output.txt";
   double **conc_mob_arr = conc[MOBILE];
   double **conc_immob_arr = conc[IMMOBILE];
   double pomoc, n;
   double el_por_m = por_m->value(ppelm->centre(), ppelm->element_accessor());
   double el_por_imm = por_imm->value(ppelm->centre(), ppelm->element_accessor());
   double el_phi = phi->value(ppelm->centre(), ppelm->element_accessor());

   vystupni_soubor = fw_chem;
   //==================================================================
   // ----------- ALOKACE POLE PRO KONCENTRACE Z FLOWA ----------------
   //==================================================================

   ASSERT(P_lat != NULL,"\nP_lat NENI ALOKOVANE\n");

  //==================================================================
   // ----------------- NEJPRVE PRO MOBILNI PORY ----------------------
   //==================================================================
   n = (el_por_m) / (1 - el_por_m); //asi S/V jako zze splocha

   switch (ppelm->dim()) { //objem se snad na nic nepouzzi:va:
	  case 1 :
	  case 2 :
	  case 3 : pomoc = ppelm->measure() *
                     cross_section->value(ppelm->centre(), ppelm->element_accessor() ) * 
                     el_por_m;
             break;
	default:
	  pomoc = 1.0;
   }

   G_prm.objem = pomoc; //objem * mobilni: porozita
   G_prm.splocha = (pomoc / el_por_m) * (el_phi) * (1 - el_por_m - el_por_imm);
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
   switch (ppelm->dim()) { //objem se snad na nic nepouzzi:va:
	  case 1 :
	  case 2 :
	  case 3 : pomoc = ppelm->measure() *
                     cross_section->value(ppelm->centre(), ppelm->element_accessor() ) * 
                     el_por_imm;
             break;
	default:
	  pomoc = 1.0;
    }
   G_prm.objem = pomoc;
   G_prm.splocha = (pomoc / el_por_imm) * (1 - el_phi) * (1 - el_por_m - el_por_imm);
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
