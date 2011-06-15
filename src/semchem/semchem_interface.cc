//---------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include "../system/system.hh"
#include "che_semchem.h"
#include "semchem_interface.hh"
#include "transport.h"

using namespace std;

//---------------------------------------------------------------------------
//  GLOBALNI PROMENNE
//---------------------------------------------------------------------------
struct TS_prm	G_prm;
struct TS_lat 	*P_lat;
struct TS_che	*P_che;

//---------------------------------------------------------------------------
Semchem_interface::Semchem_interface(void)
	:semchem_on(false)
{
  FILE *fw_chem;

  //fw_chem = fopen("vystup.txt","w"); fclose(fw_chem); //makes chemistry output file clean, before transport is computed
  fw_chem = fopen("vystup.txt","w"); fclose(fw_chem); //makes chemistry output file clean, before transport is computed
  this->semchem_on = OptGetBool("Semchem_module", "Compute_reactions", "no");
  if(semchem_on == true) ctiich();
  return;
}

//---------------------------------------------------------------------------
//                 FOR-LOOP CALLING FOR ALL ELEMENTS
//---------------------------------------------------------------------------
void Semchem_interface::compute_one_step(bool porTyp, double time_step, ElementIter ppelm, double ***conc)
{
	int nr_of_elements = 10; //this is a cheat, because I do not know how to get number of elements from this place in source code

	for (int loc_el = 0; loc_el < nr_of_elements; loc_el++)
	{
	   START_TIMER("semchem_step");
	   if(this->semchem_on == true) this->compute_reaction(porTyp, time_step, ppelm, loc_el, conc[MOBILE], conc[IMMOBILE]);
	   END_TIMER("semchem_step");
	}
}

//---------------------------------------------------------------------------
//                 FUNKCE NA VYPOCET CHEMIE PRO ELEMENT V TRANSPORTU
//---------------------------------------------------------------------------
void Semchem_interface::compute_reaction(bool porTyp, double time_step, ElementIter ppelm, int poradi, double **conc_mob_arr, double **conc_immob_arr)
{
  FILE *fw, *stream, *fr;
   int i, j; //, poradi;
   int krok;
   int poc_krok;
   double celkova_molalita;
   char vystupni_soubor[] = "vystup.txt";
   double pomoc, n;
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
	 //xprintf(Msg,"\nmolarni hmotnost: %f, latka %d, poradi %d, koncentrace z flowa %f,\n",P_lat[i].M,i,poradi,problem->transport->pconc[mobile][i][poradi]);
	 P_lat[i].m0 = (double)((conc_mob_arr[i][poradi])) / (P_lat[i].M);
	 /*if ((problem->transport->sorption == true))
	 {
		 P_lat[i].m0_sorb = (double)((sorb_mob_arr[i][poradi]) / (P_lat[i].M));
	 }
	 else
	 {
		 P_lat[i].m0_sorb = 0.0;
	 }*/
	 //xprintf(Msg,"\n %d POCATECNI KONCENTRACE transport->pconc[mobile][%d][%d]: %f, %f, sorbovane %f \n nejen G_prm.objem: %f, %f, %f\n",poradi,i,poradi,conc_mob_arr[i][poradi],P_lat[i].m0,P_lat[i].m0_sorb,ppelm->volume,pomoc,G_prm.objem);
	 if (G_prm.vypisy>1)
	 {
	  //xprintf(Msg,"\n slozka %d, mol. hmotnost %f, jednotka g/l, z transportu rslo %f, v chemii molalita %f\n", i, P_lat[i].M, problem->transport->pconc[mobile][i][poradi], P_lat[i].m0);
	 }
	 celkova_molalita += (P_lat[i].m0);
   }
   xprintf(Msg,"\nCelkova molalita: %f\n",celkova_molalita);   
   G_prm.deltaT = time_step/G_prm.cas_kroku; // dosazeni "spravneho" casoveho kroku

    //-----------------------VYPOCET CHEMIE----------------------------
   if (celkova_molalita > 1e-16)
   {
	  if (G_prm.vypisy==1)
	  {
		fw = fopen(vystupni_soubor, "a");
		fprintf(fw,"\nie = %d, dt = %lf\n", poradi, time_step);
		fclose(fw);
		che_nadpis__soubor(vystupni_soubor);
	  }
	  else if (G_prm.vypisy>1)
	  {
		fw = fopen(vystupni_soubor, "a");
		che_outpocp_soubor(fw);
		fclose(fw);
	  }
	  for (krok = 1; krok <= G_prm.cas_kroku; krok++)
	  {
		 fw = fopen(vystupni_soubor, "a");
		 if (G_prm.vypisy>1)
		 {
			fprintf(fw,"\n..............................................................\ncasovy krok c. %d, cas %f:\n", krok, krok*G_prm.deltaT);
		 }
		 else if (G_prm.vypisy==1)
		 {
			fprintf(fw,"\n%d\t%f", krok, krok*G_prm.deltaT);
		 }
		 fclose(fw);
		 che_pocitej_soubor(vystupni_soubor, &poc_krok);
		 for (i = 0; i < G_prm.pocet_latekvefazi; i++) {
		   P_lat[i].m0 = P_lat[i].m;
		   //if(problem->transport->sorption == true) P_lat[i].m0_sorb = P_lat[i].m_sorb;
		 }
		 if (G_prm.vypisy>1)
		 {
			che_vypis_soubor(vystupni_soubor);
		 }
		 else if (G_prm.vypisy==1)
		 {
			che_vypis__soubor(vystupni_soubor);
		 }
		 che_presun_poc_p_();
	  }

    //------------PREDANI VSECH INFORMACI ZE SEMCHEMU DO FLOWA---------
	  for ( i=0 ; i<G_prm.pocet_latekvefazi; i++)
	  {
			conc_mob_arr[i][poradi] = (double)(P_lat[i].m0 * P_lat[i].M);
			/*if (problem->transport->sorption == true)
			{
			   sorb_mob_arr[i][poradi] = (double)(P_lat[i].m0_sorb * (P_lat[i].M ));
			}*/
	  }
	}
   xprintf(Msg,"\n chemie je v pulce\n");

    //==================================================================
    // ----------------- POTE PRO IMOBILNI PORY ------------------------
    //==================================================================
   if (porTyp == true) {
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

		 /*if(problem->transport->sorption == true)
		 {
			 P_lat[i].m0_sorb = (sorb_immob_arr[i][poradi] / P_lat[i].M);
		 }*/
	  celkova_molalita += P_lat[i].m0;
   }

    //-----------------------VYPOCET CHEMIE----------------------------
  if (celkova_molalita > 1e-16)
   {
	  if (G_prm.vypisy==1) che_nadpis__soubor(vystupni_soubor);

	  for (krok = 1; krok<=G_prm.cas_kroku; krok++)
	  {
		 fw = fopen(vystupni_soubor, "a");
		 if (G_prm.vypisy>1)
			fprintf(fw,"\n..............................................................\ncasovy krok c. %d, cas %f:\n", krok, krok*G_prm.deltaT);
		 else if (G_prm.vypisy==1)
			fprintf(fw,"\n%d\t%f", krok, krok*G_prm.deltaT);
		 fclose(fw);
		 che_pocitej_soubor(vystupni_soubor,&poc_krok);
		 if (G_prm.vypisy>1)
			che_vypis_soubor(vystupni_soubor);
		 else if (G_prm.vypisy==1)
			che_vypis__soubor(vystupni_soubor);
		 che_presun_poc_p_();
	  }

    //------------PREDANI VSECH INFORMACI ZE SEMCHEMU DO FLOWA---------
	  for ( i=0 ; i<G_prm.pocet_latekvefazi ; i++)
	  {
			conc_immob_arr[i][poradi] = (P_lat[i].m0 * P_lat[i].M); //m nebo m0, co?
			/*if(problem->transport->sorption == true)
			{
			  sorb_immob_arr[i][poradi] == P_lat[i].m0 * P_lat[i].M;
			}*/
	  //xprintf(Msg,"\n transport->conc[immobile][i][%d]: %f \n",i,problem->transport->conc[immobile][i][poradi]);
	  }
      }
    }else{
	  xprintf(Msg,"\n chemie pro imobilni: po:ry se nepocita\n");
    }
   //} //closing bracket for the loop FOR_ELEMENTS
   xprintf(Msg,"\n skoncila chemie\n");
}

double Semchem_interface::set_timestep(double new_timestep)
{
	return new_timestep; //implementation will folow
}

