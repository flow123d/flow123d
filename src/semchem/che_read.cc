/*!
 *
﻿ * Copyright (C) 2015 Technical University of Liberec.  All rights reserved.
 * 
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License version 3 as published by the
 * Free Software Foundation. (http://www.gnu.org/licenses/gpl-3.0.en.html)
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * 
 * @file    che_read.cc
 * @brief   
 */

#include "che_semchem.h"
#include "io/read_ini.h"
#include <cstring>

extern struct TS_prm    G_prm;
extern struct TS_lat    *P_lat;
extern struct TS_che    *P_che;

// Cte CHEMIE-OBECNE ze souboru parametru .ich
void ctiich_obecne( void )
{
        int     i;
        const char *section = "Semchem_module";
        char  buffer[ 1024 ];
        char *pString;

        G_prm.pocet_latekvefazi = OptGetInt("Transport", "N_substances", NULL );
        if (G_prm.pocet_latekvefazi < 1){
          printf("\nNumber of aqueous species must be higher then 1.");
          exit(133);
        }
/*------------------------------------------------------------------*/
        G_prm.T = OptGetDbl(section,"Temperature","0.0");
        if ( G_prm.T <= 0.0 )
        {
        printf("\nTeplota musi byt kladna!");
      exit(133);
        }
/*------------------------------------------------------------------*/
        G_prm.TGf = OptGetDbl(section,"Temperature_Gf","-1.0");
        if ( G_prm.TGf == -1.0 )
        {
      G_prm.TGf = G_prm.T;
        }
        if ( G_prm.T <= 0.0 )
        {
        printf("\nTeplota pro zadani dGf musi byt kladna!");
      exit(133);
        }
/*------------------------------------------------------------------*/
        G_prm.epsilon = OptGetDbl(section,"Epsilon","0.0");
        if ( G_prm.epsilon<= 0.0 )
        {
        printf("\nEpsilon musi byt kladne!");
      exit(133);
        }
/*------------------------------------------------------------------*/
   pString = strcpy(buffer,OptGetStr(section, "Error_norm_type", "Absolute"));
   if ( pString == NULL )
   {
      printf("\nChybi typ normy!");
      exit(133);
   }
   G_prm.abs_norma = -1;
   if ( strcmp( buffer, "relative" ) == 0 ) G_prm.abs_norma = 0;
   if ( strcmp( buffer, "Relative" ) == 0 ) G_prm.abs_norma = 0;
   if ( strcmp( buffer, "rel" ) == 0 ) G_prm.abs_norma = 0;
   if ( strcmp( buffer, "Rel" ) == 0 ) G_prm.abs_norma = 0;
   if ( strcmp( buffer, "r" ) == 0 ) G_prm.abs_norma = 0;
   if ( strcmp( buffer, "R" ) == 0 ) G_prm.abs_norma = 0;
   if ( strcmp( buffer, "0" ) == 0 ) G_prm.abs_norma = 0;
   if ( strcmp( buffer, "absolute" ) == 0 ) G_prm.abs_norma = 1;
   if ( strcmp( buffer, "Absolute" ) == 0 ) G_prm.abs_norma = 1;
   if ( strcmp( buffer, "abs" ) == 0 ) G_prm.abs_norma = 1;
   if ( strcmp( buffer, "Abs" ) == 0 ) G_prm.abs_norma = 1;
   if ( strcmp( buffer, "a" ) == 0 ) G_prm.abs_norma = 1;
   if ( strcmp( buffer, "A" ) == 0 ) G_prm.abs_norma = 1;
   if ( strcmp( buffer, "1" ) == 0 ) G_prm.abs_norma = 1;
   if (G_prm.abs_norma == -1)
   {
      printf("\nTyp normy neni platny!");
      exit(133);
   }
   G_prm.omega = 1.0;
/*------------------------------------------------------------------*/
   pString = strcpy(buffer,OptGetStr(section,"Scaling","No"));
   if ( pString == NULL )
   {
      printf("\nChybi definice skalovani!");
      exit(133);
   }
   G_prm.skaluj_matici = -1;
   if ( strcmp( buffer, "N" ) == 0 ) G_prm.skaluj_matici = 0;
   if ( strcmp( buffer, "Ne" ) == 0 ) G_prm.skaluj_matici = 0;
   if ( strcmp( buffer, "No" ) == 0 ) G_prm.skaluj_matici = 0;
   if ( strcmp( buffer, "0" ) == 0 ) G_prm.skaluj_matici = 0;
   if ( strcmp( buffer, "A" ) == 0 ) G_prm.skaluj_matici = 1;
   if ( strcmp( buffer, "Ano" ) == 0 ) G_prm.skaluj_matici = 1;
   if ( strcmp( buffer, "Y" ) == 0 ) G_prm.skaluj_matici = 1;
   if ( strcmp( buffer, "Yes" ) == 0 ) G_prm.skaluj_matici = 1;
   if ( strcmp( buffer, "1" ) == 0 ) G_prm.skaluj_matici = 1;
   if (G_prm.skaluj_matici == -1)
   {
      printf("\nSkalovani matice neni platne!");
      exit(133);
   }
/*-----------------------------------------------*/
        G_prm.Afi = OptGetDbl(section,"Param_Afi","-1.0");
        if ( G_prm.Afi < 0.0 )
        {
        printf("\nAfi musi byt nezaporne!");
      exit(133);
   }
/*------------------------------------------------------------------*/
        G_prm.b = OptGetDbl(section,"Param_b","0.0");
        if ( G_prm.b <= 0.0 )
        {
        printf("\nb musi byt kladne!");
      exit(133);
   }
/*------------------------------------------------------------------*/
        G_prm.cas_kroku = OptGetInt(section,"Time_steps","0");
        if ( G_prm.cas_kroku <= 0 )
        {
        printf("\nPocet casovych kroku musi byt kladny!");
      exit(133);
   }
/*------------------------------------------------------------------*/
        G_prm.vypisy = OptGetInt(section,"Output_precission","0");
/*------------------------------------------------------------------*/
        i = OptGetInt(section, "Number_of_further_species","0");
        if ( i < 0 )
        {
        printf("\nPocet dalsich latek nesmi byt zaporny!");
      exit(133);
   }
   G_prm.pocet_latek = i+G_prm.pocet_latekvefazi;
   if (G_prm.pocet_latek>MAX_POC_LATEK)
   {
           printf ("Celkovy pocet latek muze byt maximalne %d!\n", MAX_POC_LATEK);
           exit(121);
   }
/*------------------------------------------------------------------*/
        G_prm.deleni_RK = OptGetInt(section,"Slow_kinetics_substeps","1");
        if ( G_prm.deleni_RK < 1 )
        {
        printf("\nPocet kroku pomale kinetiky musi byt kladny!");
      exit(133);
   }
}
/********************************************************************/
/*                      Cte LATKY_VE_FAZI ze souboru parametru .ich */
/********************************************************************/
void ctiich_latkyvefazi( void )
{
   char  buffer[1024];
   int   j;
   char  nazev[30], *pom_buf;
   const char *separators = " ,\t";

// Alokace seznamu latek
   P_lat = (TS_lat *)malloc( (G_prm.pocet_latek)*sizeof( TS_lat ) );
   if ( P_lat == NULL )
   {
           printf ("Malo pameti!\n");
           exit(0);
   }
// Nacteni obsahu seznamu latek
/*-----------------------------------------------*/
   sprintf( nazev, "Aqueous_species" );
    {
       for (j=0; j<G_prm.pocet_latekvefazi; j++)
       {
          strcpy(P_lat[j].nazev,"");
       }
    }
/*-----------------------------------------------*/
   strcpy(buffer,OptGetStr(nazev,"dGf","<NeplatnyNazev>"));
   pom_buf = strtok( buffer, separators );
   if ( strcmp( buffer, "<NeplatnyNazev>" ) == 0 )
   {
      pom_buf = NULL;
   }
   for (j=0; j<G_prm.pocet_latekvefazi; j++)
   {
      if ( pom_buf == NULL )
      {
         printf("\nChybi dGf %d. latky ve fazi!", j+1);
         exit(133);
      }
      P_lat[j].dGf = atof(pom_buf);
        printf("\n P_lat[%d].dGf %f\n",j,P_lat[j].dGf);
      pom_buf = strtok( NULL, separators );
   }
   if ( pom_buf != NULL )
   {
      printf("\nPrilis mnoho dGf pro latky ve fazi!");
      exit(133);
   }
/*-----------------------------------------------*/
   strcpy(buffer,OptGetStr(nazev,"dHf","<NeplatnyNazev>"));
   pom_buf = strtok( buffer, separators );
   if ( strcmp( buffer, "<NeplatnyNazev>" ) == 0 )
   {
      for (j=0; j<G_prm.pocet_latekvefazi; j++)
      {
         P_lat[j].dHf = 0.0;
      }
      pom_buf = NULL;
   }
   else for (j=0; j<G_prm.pocet_latekvefazi; j++)
   {
      if ( pom_buf == NULL )
      {
         printf("\nChybi dHf %d. latky ve fazi!", j+1);
         exit(133);
      }
      P_lat[j].dHf = atof(pom_buf);
      pom_buf = strtok( NULL, separators );
   }
   if ( pom_buf != NULL )
   {
      printf("\nPrilis mnoho dHf pro latky ve fazi!");
      exit(133);
   }
/*-----------------------------------------------*/
   strcpy(buffer,OptGetStr(nazev,"Molar_mass","<NeplatnyNazev>"));
   pom_buf = strtok( buffer, separators );
   if ( strcmp( buffer, "<NeplatnyNazev>" ) == 0 )
   {
      pom_buf = NULL;
   }
      for (j=0; j<G_prm.pocet_latekvefazi; j++)
      {
         if ( pom_buf == NULL )
         {
            printf("\nChybi molarni hmotnost %d. latky ve fazi!", j+1);
            exit(133);
         }
                 P_lat[j].M = atof(pom_buf);
         pom_buf = strtok( NULL, separators );
      }
      if ( pom_buf != NULL )
      {
         printf("\nPrilis mnoho molarnich hmotnosti pro latky ve fazi!");
         exit(133);
      }
/*-----------------------------------------------*/
   strcpy(buffer,OptGetStr(nazev,"El_charge","<NeplatnyNazev>"));
   pom_buf = strtok( buffer, separators );
   if ( strcmp( buffer, "<NeplatnyNazev>" ) == 0 )
   {
      pom_buf = NULL;
   }
   for (j=0; j<G_prm.pocet_latekvefazi; j++)
   {
      if ( pom_buf == NULL )
      {
         printf("\nChybi naboj %d. latky ve fazi!", j+1);
         exit(133);
      }
      P_lat[j].Q = atoi(pom_buf);
      pom_buf = strtok( NULL, separators );
   }
   if ( pom_buf != NULL )
   {
      printf("\nPrilis mnoho naboju pro latky ve fazi!");
      exit(133);
   }
/*-----------------------------------------------*/
          for (j=0; j<G_prm.pocet_latekvefazi; j++)
          {
                 P_lat[j].typ_sorpce = 0; //Sorption has been suppressed in semchem module.
          }
/*-----------------------------------------------*/
}
// Cte DALSI_LATKY ze souboru parametru .ini
void ctiich_dalsilatky( void )
{
   char  buffer[1024];
   int   j;
   char  nazev[30], *pom_buf;
   const char* separators = " ,\t";

   if (G_prm.pocet_latekvefazi == 0)
   {
      return;
   }
// Nacteni obsahu seznamu latek
   sprintf( nazev, "Further_species" );
   strcpy(buffer,OptGetStr(nazev,"Specie_name","<NeplatnyNazev>"));
   pom_buf = strtok( buffer, separators );
   if ( strcmp( buffer, "<NeplatnyNazev>" ) == 0 )
   {
      for (j=G_prm.pocet_latekvefazi;j<G_prm.pocet_latek; j++)
      {
         sprintf( P_lat[j].nazev, "Dalsi_latka_%d", j+1-G_prm.pocet_latek );
      }
   }
   else
   {
      for (j=G_prm.pocet_latekvefazi;j<G_prm.pocet_latek; j++)
      {
         if ( pom_buf == NULL )
         {
            printf("\nChybi nazev %d. dalsi latky!", j+1-G_prm.pocet_latek);
            exit(133);
         }
         strcpy (P_lat[j].nazev, pom_buf);
         pom_buf = strtok( NULL, separators );
      }
      if ( pom_buf != NULL )
      {
         printf("\nPrilis mnoho nazvu dalsich latek!");
         exit(133);
      }
   }
/*-----------------------------------------------*/
   strcpy(buffer,OptGetStr(nazev,"dGf","<NeplatnyNazev>"));
   pom_buf = strtok( buffer, separators );
   if ( strcmp( buffer, "<NeplatnyNazev>" ) == 0 )
   {
      pom_buf = NULL;
   }
   for (j=G_prm.pocet_latekvefazi;j<G_prm.pocet_latek; j++)
   {
      if ( pom_buf == NULL )
      {
         printf("\nChybi dGf %d. dalsi latky!", j+1-G_prm.pocet_latek);
         exit(133);
      }
      P_lat[j].dGf = atof(pom_buf);
      pom_buf = strtok( NULL, separators );
   }
   if ( pom_buf != NULL )
   {
      printf("\nPrilis mnoho dGf pro dalsi latky!");
      exit(133);
   }
/*-----------------------------------------------*/
   strcpy(buffer,OptGetStr(nazev,"dHf","<NeplatnyNazev>"));
   pom_buf = strtok( buffer, separators );
   if ( strcmp( buffer, "<NeplatnyNazev>" ) == 0 )
   {
      for (j=G_prm.pocet_latekvefazi;j<G_prm.pocet_latek; j++)
      {
         P_lat[j].dHf = 0.0;
      }
      pom_buf = NULL;
   }
   else for (j=G_prm.pocet_latekvefazi;j<G_prm.pocet_latek; j++)
   {
      if ( pom_buf == NULL )
      {
         printf("\nChybi dHf %d. dalsi latky!", j+1-G_prm.pocet_latek);
         exit(133);
      }
      P_lat[j].dHf = atof(pom_buf);
      pom_buf = strtok( NULL, separators );
   }
   if ( pom_buf != NULL )
   {
      printf("\nPrilis mnoho dHf pro dalsi latky!");
      exit(133);
   }
/*-----------------------------------------------*/
   strcpy(buffer,OptGetStr(nazev,"Molar_mass","<NeplatnyNazev>"));
   pom_buf = strtok( buffer, separators );
   if ( strcmp( buffer, "<NeplatnyNazev>" ) == 0 )
   {
      pom_buf = NULL;
   }
   for (j=G_prm.pocet_latekvefazi;j<G_prm.pocet_latek; j++)
   {
      if ( pom_buf == NULL )
      {
         printf("\nChybi molarni hmotnost %d. dalsi latky!", j+1-G_prm.pocet_latek);
         exit(133);
      }
      P_lat[j].M = atof(pom_buf);
      pom_buf = strtok( NULL, separators );
   }
   if ( pom_buf != NULL )
   {
      printf("\nPrilis mnoho molarnich hmotnosti pro dalsi latky!");
      exit(133);
   }
/*-----------------------------------------------*/
   strcpy(buffer,OptGetStr(nazev,"Activity","<NeplatnyNazev>"));
   pom_buf = strtok( buffer, separators );
   if ( strcmp( buffer, "<NeplatnyNazev>" ) == 0 )
   {
      pom_buf = NULL;
   }
   for (j=G_prm.pocet_latekvefazi;j<G_prm.pocet_latek; j++)
   {
      if ( pom_buf == NULL )
      {
         printf("\nChybi aktivita %d. dalsi latky!", j+1-G_prm.pocet_latek);
         exit(133);
      }
      P_lat[j].aktivita = atof(pom_buf);
      pom_buf = strtok( NULL, separators );
   }
   if ( pom_buf != NULL )
   {
      printf("\nPrilis mnoho aktivit pro dalsi latky!");
      exit(133);
   }
}

// Cte REAKCE ze souboru parametru .ich
void ctiich_reakce( void )
{
   char  buffer[1024];
   int   i,j;
   char  nazev[30], *pom_buf;
   const char* separators = " ,\t";
   TS_che pom_che;
   double max_poloc_rozp;

// Zjisteni delky seznamu reakci
   for ( i = 1;; ++i )
   {
      sprintf( nazev, "Reaction_%d", i );
      strcpy(buffer,OptGetStr(nazev,"Reaction_type","<NeplatnyNazev>")); printf("probehlo zjisteni typu reakce %d: %s\n",i, buffer);
      if ( strcmp( buffer, "<NeplatnyNazev>" ) == 0 ) break;
   }
   G_prm.celkovy_pocet_reakci = i-1;
        if ( i==1 )
        {
        printf("\nNeni definovana zadna reakce!");
      exit(133);
   }
// Alokace seznamu reakci
   P_che = (TS_che *)malloc( (G_prm.celkovy_pocet_reakci)*sizeof( TS_che ) );
   if ( P_che == NULL )
   {
           printf ("Malo pameti!\n");
           exit(0);
   }
// Nacteni obsahu seznamu reakci
        G_prm.pocet_reakci_pro_matici = 0;
   G_prm.pocet_rozpadu = 0;
   G_prm.pocet_pom_kin = 0;
   max_poloc_rozp = 0.0;
   for (i=0; i<G_prm.celkovy_pocet_reakci; i++)
   {
      sprintf( nazev, "Reaction_%d", i+1 );
      sprintf(P_che[i].nazev,"Reaction_%d",i+1);
/*-----------------------------------------------*/
      strcpy(buffer,OptGetStr(nazev,"Reaction_type","<NeplatnyTyp>"));
      pom_buf = strtok( buffer, separators );
      P_che[i].typ_reakce = -888;
      if ( strcmp( pom_buf, "Radioactive_decay" ) == 0 ) P_che[i].typ_reakce = 4;
      if ( strcmp( pom_buf, "radioactive_decay" ) == 0 ) P_che[i].typ_reakce = 4;
      if ( strcmp( pom_buf, "Decay" ) == 0 ) P_che[i].typ_reakce = 4;
      if ( strcmp( pom_buf, "decay" ) == 0 ) P_che[i].typ_reakce = 4;
      if ( strcmp( pom_buf, "RD" ) == 0 ) P_che[i].typ_reakce = 4;
      if ( strcmp( pom_buf, "rd" ) == 0 ) P_che[i].typ_reakce = 4;
      if ( strcmp( pom_buf, "4" ) == 0 ) P_che[i].typ_reakce = 4;
      if ( strcmp( pom_buf, "Slow_kinetics" ) == 0 ) P_che[i].typ_reakce = 3;
      if ( strcmp( pom_buf, "slow_kinetics" ) == 0 ) P_che[i].typ_reakce = 3;
      if ( strcmp( pom_buf, "3" ) == 0 ) P_che[i].typ_reakce = 3;
      if ( strcmp( pom_buf, "Kinetics" ) == 0 ) P_che[i].typ_reakce = 1;
      if ( strcmp( pom_buf, "kinetics" ) == 0 ) P_che[i].typ_reakce = 1;
      if ( strcmp( pom_buf, "1" ) == 0 ) P_che[i].typ_reakce = 1;
      if ( strcmp( pom_buf, "Equilibrium" ) == 0 ) P_che[i].typ_reakce = 0;
      if ( strcmp( pom_buf, "equilibrium" ) == 0 ) P_che[i].typ_reakce = 0;
      if ( strcmp( pom_buf, "0" ) == 0 ) P_che[i].typ_reakce = 0;
      if (P_che[i].typ_reakce == -888)
      {
         printf("\nThe type of reaction nr. %d is not valid!", i);
         exit(133);
      }
      if (P_che[i].typ_reakce==0)
      {
         G_prm.pocet_rovnovah++;
      }
      if (P_che[i].typ_reakce==1)
      {
         G_prm.pocet_kinetik++;
      }
      if (P_che[i].typ_reakce==3)
      {
         G_prm.pocet_pom_kin++;
      }
      if (P_che[i].typ_reakce==4)
      {
         G_prm.pocet_rozpadu++;
      }
      G_prm.pocet_reakci_pro_matici = G_prm.pocet_rovnovah+G_prm.pocet_kinetik;
/*------------------------------------------------------------------*/
      strcpy(buffer,OptGetStr(nazev,"Stoichiometry",""));
      pom_buf = strtok( buffer, separators );
      for (j = 0; j<G_prm.pocet_latek; j++)
      {
         if ( pom_buf == NULL )
         {
            printf("\nChybi %d. stechometricky koeficient v %d. rovnici!", j+1, i+1);
            exit(133);
         }
                 P_che[i].stech_koef_p[j] = atoi( pom_buf );
                 printf("\nP_che[%d].stech_koef_p[%d]: %d",i,j,P_che[i].stech_koef_p[j]);
         pom_buf = strtok( NULL, separators );
      }
      if ( pom_buf != NULL )
      {
         printf("\nV %d. rovnici je prilis mnoho stechiometrickych koeficientu!", i+1);
         exit(133);
      }
/*------------------------------------------------------------------*/
      if ((P_che[i].typ_reakce==1)||(P_che[i].typ_reakce==3))
      {
         P_che[i].K = OptGetDbl(nazev,"Kinetic_constant","0.0");
         if(P_che[i].typ_reakce==1)printf("\nKineticka konstanta v %d. rovnici ma hodnotu %f", i+1, P_che[i].K);
         if ( P_che[i].K <= 0.0 )
         {
            printf("\nKineticka konstanta v %d. rovnici neni kladna!", i+1);
            exit(133);
         }
/*------------------------------------------------------------------*/
         strcpy(buffer,OptGetStr(nazev,"Order_of_reaction",""));
         pom_buf = strtok( buffer, separators );
         for (j = 0; j<G_prm.pocet_latek; j++)
         {
            if(j < G_prm.pocet_latekvefazi)
            {
                if( pom_buf == NULL )
                {
               printf("\nChybi %d. mocnina pro kinetiku %d. rovnice!", j+1, i+1);
               exit(133);
                }
            P_che[i].exponent[j] = atof( pom_buf );
            printf("\nP_che[%d].exponent[%d]: %f",i,j,P_che[i].exponent[j]);
            pom_buf = strtok( NULL, separators );
           }
         }
         if ( pom_buf != NULL )
         {
            printf("\nV %d. rovnici je prilis mnoho exponentu pro kinetiku!", i+1);
            exit(133);
         }
      }
/*------------------------------------------------------------------*/
      if (P_che[i].typ_reakce==0)
      {
         P_che[i].K = OptGetDbl(nazev, "Equilibrium_constant","-1.0");
      }
/*------------------------------------------------------------------*/
   }
        if (G_prm.celkovy_pocet_reakci>1)
   {
      for (j=0; j<G_prm.celkovy_pocet_reakci; j++)
      {
         for (i=0; i<G_prm.celkovy_pocet_reakci-1; i++)
         {
            if (che_poradi(P_che[i].typ_reakce,max_poloc_rozp,P_che[i].K) > che_poradi(P_che[i+1].typ_reakce,max_poloc_rozp,P_che[i+1].K))
            {
                pom_che = P_che[i];
               P_che[i] = P_che[i+1];
               P_che[i+1] = pom_che;
            }
         }
      }
   }
}

// Cte soubor parametru chemie .ICH
void ctiich( void )
{
        printf("Cteni parametru (%s): \n", G_prm.jmeno_ich );
        ctiich_obecne(); printf("probehla funkce ctiich_obecne()\n");
        ctiich_latkyvefazi(); printf("probehla funkce ctiich_latkyvefazi()\n");
        ctiich_dalsilatky(); printf("probehla funkce ctiich_dalsilatky()\n");
        ctiich_reakce(); printf("probehla funkce ctiich_reakce()\n");
        printf("O.K.\n" );
}
