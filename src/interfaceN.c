//---------------------------------------------------------------------------

#pragma hdrstop
#include "interfaceN.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <petscvec.h>
#include "materials.hh"
#include "transport.h"
#include "mesh.h"
#include "problem.h"
#include "system.hh"
#include "mesh/ini_constants_mesh.hh"
#include "constantdb.h"

//--------------pro semchem--------------------------------------------------
// void che_nadpis__soubor(char *soubor);
// void che_outpocp_soubor(FILE *fw);
// void che_pocitej_soubor(char *soubor, int *poc_krok);
// void che_vypis_soubor(char *soubor);
// void che_presun_poc_p_(void);
// void che_vypis__soubor(char *soubor);

//---------------------------------------------------------------------------
//  GLOBALNI PROMENNE
//---------------------------------------------------------------------------
struct TS_prm	G_prm;
struct TS_lat 	*P_lat;
struct TS_che	*P_che;

//---------------------------------------------------------------------------
void priprav(char *jmeno){
  strcpy(G_prm.jmeno_ich, jmeno);
  ctiich();
  return;
}

//---------------------------------------------------------------------------
void kontrola(){
   FILE *fw;
   int i,j;
 
   fw = fopen("kontrola.txt","w");
   fprintf(fw,"\nKONTROLA\n-------");
   fprintf(fw,"\n%s",G_prm.jmeno_ich);
   fprintf(fw,"\nG_prm.pocet_latek: %d",G_prm.pocet_latek);
   fprintf(fw,"\nG_prm.pocet_latekvefazi: %d",G_prm.pocet_latekvefazi);
   fprintf(fw,"\nG_prm.celkovy_pocet_reakci: %d",G_prm.celkovy_pocet_reakci);
   fprintf(fw,"\nG_prm.pocet_reakci_pro_matici: %d",G_prm.pocet_reakci_pro_matici);
   fprintf(fw,"\nG_prm.pocet_rovnovah: %d",G_prm.pocet_rovnovah);
   fprintf(fw,"\nG_prm.pocet_kinetik: %d",G_prm.pocet_kinetik);
   fprintf(fw,"\nG_prm.pocet_pom_kin: %d",G_prm.pocet_pom_kin);
   fprintf(fw,"\nG_prm.pocet_rozpadu: %d",G_prm.pocet_rozpadu);
   fprintf(fw,"\nG_prm.T: %f",G_prm.T);
   fprintf(fw,"\nG_prm.TGf: %f",G_prm.TGf);
   fprintf(fw,"\nG_prm.Afi: %f",G_prm.Afi);
   fprintf(fw,"\nG_prm.b: %f",G_prm.b);
   fprintf(fw,"\nG_prm.epsilon: %f",G_prm.epsilon);
   fprintf(fw,"\nG_prm.omega: %f",G_prm.omega);
   fprintf(fw,"\nG_prm.deltaT: %f",G_prm.deltaT);
   fprintf(fw,"\nG_prm.objem: %f",G_prm.objem);
   fprintf(fw,"\nG_prm.splocha: %f",G_prm.splocha);
   fprintf(fw,"\nG_prm.vypisy: %d",G_prm.vypisy);
   for (j = 0; j < G_prm.celkovy_pocet_reakci; j++){
 	 for(i=0;i<G_prm.pocet_latek;i++){
 	  fprintf(fw,"\n %s, P_che[%d].stech_koef_p[%d]: %d",P_che[j].nazev,j,i,P_che[j].stech_koef_p[i]);
 	  //if(i < (G_prm.pocet_latek - 1)) 
//  	    xprintf(Msg,", P_che[%d].exponent[%d]: %d\n mola:rni: hmotnost: %f",j,i,P_che[j].exponent[i],P_lat[i].M);
 	  //else 
 	    fprintf(fw,"\n");
 	}
   }
//    fprintf(fw,"parametr sorpce %f",P_lat[0].param_sorpce[0]);
   fclose(fw);
   return;
}

//---------------------------------------------------------------------------
//                 FUNKCE NA VYPOCET CHEMIE PRO ELEMENT V TRANSPORTU
//---------------------------------------------------------------------------
//void che_vypocetchemie(struct Problem *problem, Vec *conc_mob_arr, Vec *conc_immob_arr, Vec *sorb_mob_arr, Vec *sorb_immob_arr)
//void che_vypocetchemie(struct Problem *problem, PetscScalar **conc_mob_arr, PetscScalar **conc_immob_arr, PetscScalar **sorb_mob_arr, PetscScalar **sorb_immob_arr)
void che_vypocetchemie(struct Problem *problem, double **conc_mob_arr, double **conc_immob_arr, double **sorb_mob_arr, double **sorb_immob_arr)
{
   Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);
  
   ElementIter ppelm;
   FILE *fw, *stream, *fr;
   int i, j, poradi; //, elm_id;
   int krok;
   int poc_krok;
   int n_elements = mesh->n_elements();
   double celkova_molalita;
   char vystupni_soubor[] = "vystup.txt";
   double pomoc, n;
   bool porTyp; 
   //double ***concentracion = problem->transport->pconc;
   //double **concentracion;
   //==================================================================
   // ----------- ALOKACE POLE PRO KONCENTRACE Z FLOWA ----------------
   //==================================================================   
   //concentracion = (double ***) malloc(sizeof(double *) * problem->transport->n_substances);
   //místo pro předání pointerů na PETSc-í vektory
   //==================================================================

   ASSERT(P_lat != NULL,"\nP_lat NENI ALOKOVANE\n");

   porTyp = problem->transport->dual_porosity;  

   //for(poradi = 0; poradi < (n_elements); poradi++){ //bezi pres celou funkci
   FOR_ELEMENTS(ppelm){
   poradi = ELEMENT_FULL_ITER(ppelm) - mesh->element.begin();
   //elm_id = mesh->epos_id[poradi]; 
   //ppelm = &mesh->element[elm_id];
   xprintf(Msg,"\n!!!! TRANSPORT %d!!!!\n",problem->transport);
   xprintf(Msg,"\nmolarni hmotnost: %f, latka %d, poradi %d\n",P_lat[0].M,0,poradi);

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
	 xprintf(Msg,"\nmolarni hmotnost: %f, latka %d, poradi %d, koncentrace z flowa %f,\n",P_lat[i].M,i,poradi,problem->transport->pconc[mobile][i][poradi]);
	 P_lat[i].m0 = (double)((conc_mob_arr[i][poradi])) / (P_lat[i].M);
	 if ((problem->transport->sorption == true))
	 {
		 P_lat[i].m0_sorb = (double)((sorb_mob_arr[i][poradi]) / (P_lat[i].M));
	 }
	 else
	 {
		 P_lat[i].m0_sorb = 0.0;
	 }
	 xprintf(Msg,"\n %d POCATECNI KONCENTRACE transport->pconc[mobile][%d][%d]: %f, %f, sorbovane %f \n nejen G_prm.objem: %f, %f, %f\n",poradi,i,poradi,conc_mob_arr[i][poradi],P_lat[i].m0,P_lat[i].m0_sorb,ppelm->volume,pomoc,G_prm.objem);
	 if (G_prm.vypisy>1)
	 {
	  xprintf(Msg,"\n slozka %d, mol. hmotnost %f, jednotka g/l, z transportu rslo %f, v chemii molalita %f\n", i, P_lat[i].M, problem->transport->pconc[mobile][i][poradi], P_lat[i].m0);
	 }
	 celkova_molalita += (P_lat[i].m0);
   }
   xprintf(Msg,"\nCelkova molalita: %f\n",celkova_molalita);   
   G_prm.deltaT = problem->transport->time_step / G_prm.cas_kroku; // dosazeni "spravneho" casoveho kroku

    //-----------------------VYPOCET CHEMIE----------------------------
   if (celkova_molalita > 1e-16)
   {
	  if (G_prm.vypisy==1)
	  {
		fw = fopen(vystupni_soubor, "a");
		fprintf(fw,"\nie = %d, dt = %lf\n", poradi, problem->transport->time_step);
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
		   if(problem->transport->sorption == true) P_lat[i].m0_sorb = P_lat[i].m_sorb;
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
			if (problem->transport->sorption == true)
			{
			   sorb_mob_arr[i][poradi] = (double)(P_lat[i].m0_sorb * (P_lat[i].M ));
			}
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

		 if(problem->transport->sorption == true)
		 {
			 P_lat[i].m0_sorb = (sorb_immob_arr[i][poradi] / P_lat[i].M);
		 }
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
			if(problem->transport->sorption == true)
			{
			  sorb_immob_arr[i][poradi] == P_lat[i].m0 * P_lat[i].M;
			}
	  xprintf(Msg,"\n transport->conc[immobile][i][%d]: %f \n",i,problem->transport->conc[immobile][i][poradi]);
	  }
      }
    }else{
	  xprintf(Msg,"\n chemie pro imobilni: po:ry se nepocita\n");
    }
   } //closing bracket for the loop FOR_ELEMENTS
   //==================================================================
   // ----------- ALOKACE POLE PRO KONCENTRACE Z FLOWA ----------------
   //==================================================================   
   // NEEEEEEE, free(concentracion); BY UVOLNILO PAMĚŤ S PETSc-ÍMI VEKTORY A TO MY NECHCEM
   //concentration = NULL; //pro jistotu
   //==================================================================
   xprintf(Msg,"\n skoncila chemie\n");
}

void kontrola2()
{
//pro pripad, ze je to nekde volany
}
#pragma package(smart_init)
