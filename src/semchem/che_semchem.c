//#define MAIN
#include "semchem/big-head2.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "system.hh"
#define _R 0.008314
#define VERZE "30-09-2005"
#define vystupni_soubor "che_out.txt"
#define PUL_POCTU_VNORENI 8
#define MAX_POC_VNEJ_CYK 1500
#define MAX_STAGNACE 400
#define MAX_POC_VNITR_CYK 2000
#define DROBNY_POSUN 1.0e-10


void che_vypis_soubor(char *soubor)
{
   int i = 0;
   FILE *fw;

   fw = fopen(soubor, "a");
    for (i=0; i<G_prm.pocet_latekvefazi; i++)
//  	   fprintf (fw, "\nmolalita rozpustene %s : %f", P_lat[i].nazev, P_lat[i].m);
 	   fprintf (fw, "\nmolalita rozpustene %d. latky: %f", i, P_lat[i].m);
    for (i=0; i<G_prm.pocet_latekvefazi; i++)
//  	   fprintf (fw, "\nmolalita sorbovane  %s : %f", P_lat[i].nazev, P_lat[i].m_sorb);
 	   fprintf (fw, "\nmolalita sorbovane  %d. latky: %f", i, P_lat[i].m_sorb);
   fclose(fw);
}

void che_vypis__soubor(char *soubor)
{
   int i = 0;
   FILE *fw;

   fw = fopen(soubor, "a");
   for (i=0; i<G_prm.pocet_latekvefazi; i++)
	   fprintf (fw,"\t%f", P_lat[i].m);
   for (i=0; i<G_prm.pocet_latekvefazi; i++)
	   fprintf (fw,"\t%f", P_lat[i].m_sorb);
   fprintf(fw,"\t%f",G_prm.objem);
   fclose(fw);
}

void che_outpocp_soubor(FILE *fw)
{
   int i = 0;

 	fprintf(fw,"\n..............................................................");
    for (i=0; i<G_prm.pocet_latekvefazi; i++)
//  	   fprintf (fw, "\npocatecni molalita rozpustene %s : %f", P_lat[i].nazev, P_lat[i].m0);
 	   fprintf (fw, "\npocatecni molalita rozpustene %d. latky: %f", i, P_lat[i].m0);
    for (i=0; i<G_prm.pocet_latekvefazi; i++)
 	   fprintf (fw, "\npocatecni molalita sorbovane  %d. latky: %f", i, P_lat[i].m0_sorb);
}

void che_outpocp__soubor(FILE *fw)
{
   int i = 0;

   for (i=0; i<G_prm.pocet_latekvefazi; i++)
	   fprintf (fw,"\t%f", P_lat[i].m0);
   for (i=0; i<G_prm.pocet_latekvefazi; i++)
	   fprintf (fw,"\t%f", P_lat[i].m0_sorb);
   fprintf(fw,"\t%f",G_prm.objem);
}

int che_Gauss ( double *matice, double *prstrana, int *hprvky, int rozmer )
{
   int i = 0;
   int j = 0;
   int k = 0;
   int prvek = 0;
   double velprvku = 0.0;

	if (G_prm.vypisy>4) printf("\nChe_Gauss");//xprintf(Msg,"\nche_Gauss: ");
   if (rozmer <1)
      return -1;
   if (rozmer ==1)
   {
      hprvky[0]=0;
      if (matice[0] == 0.0)
      {
      	printf("prusvih Newtona: nulova derivace, prava strana = %f !\n", prstrana[0]);//xprintf(Msg,"prusvih Newtona: nulova derivace, prava strana = %f !\n", prstrana[0]);
         if ( prstrana[0] == 0.0)
         {
         	prstrana[0] = 1.0e-10;						// cislo vycucane z prstu
         	return 0;
         }
         else
         {
         	prstrana[0]*= 1.0e2;	    					// cislo vycucane z prstu
            return 0;
         }
      }
      else
      {
      	prstrana[0]/=matice[0];
      }
      return 0;
   }

   for ( i=0; i<rozmer; ++i )
   {
// Vyber hlavniho prvku
      velprvku = 0.0;
      prvek = -1;
      k=0;                                        // kvuli vyberu hlavniho prvku
      for ( j=0; j<rozmer; ++j )
      {
         if (i>0)
         {
            for ( k=0; k<i; ++k )
            {
               if (hprvky[k]==j) break;
            }
         }
         if ((k==i) && (fabs(matice[j*rozmer+i])>fabs(velprvku)))
         {
            prvek = j;
            velprvku = matice[j*rozmer+i];
         }
      }
// hlavni prvek vybran v prvek
      if (prvek == -1)
      {
         printf("Prusvih v %d. sloupci\n", i);//xprintf(Msg,"Prusvih v %d. sloupci\n", i);
//         for (j=0;j<=i;++j) printf("hprvky[%ld]= %ld\n",j,hprvky[j]);
//         return -1;
         for ( j=0; j<rozmer; ++j )
         {
            if (i>0)
            {
               for ( k=0; k<i; ++k )
               {
                  if (hprvky[k]==j) break;
               }
            }
            if (k==i)
            {
					matice[j*rozmer+i]=1.0e-9;						// cislo vycucane z prstu
               prvek = j;
               velprvku = matice[j*rozmer+i];
               break;
            }
         }
      }
      hprvky[i] = prvek;
// deleni celeho radku
      matice[prvek*rozmer+i] = 1.0;
      for ( j=i+1; j<rozmer; ++j ) matice[prvek*rozmer+j]/=velprvku;
      prstrana[prvek]/=velprvku;
// odecitani od ostatnich radku
       for ( j=0; j<rozmer; ++j )
       {
         for ( k=0; k<=i; ++k )
         {
            if (hprvky[k]==j) break;
         }
         if (k>i)
         {
              for ( k=i+1; k<rozmer; ++k ) matice[j*rozmer+k]-=matice[j*rozmer+i]*matice[prvek*rozmer+k];
            prstrana[j]-=matice[j*rozmer+i]*prstrana[prvek];
            matice[j*rozmer+i]=0.0;
         }
      }
   }
// zpetny chod
   for ( i=rozmer-1; i>=0; --i )
   {
      for ( j=0; j<rozmer; ++j )
      {
         for ( k=i; k<rozmer; ++k )
         {
            if (hprvky[k]==j) break;
         }
         if (k==rozmer)
         {
            prstrana[j]-=matice[j*rozmer+i]*prstrana[hprvky[i]];
            matice[j*rozmer+i]=0.0;
         }
      }
   }
	if (G_prm.vypisy>4) printf("o.k. (che_Gauss)");//xprintf(Msg,"o.k. (che_Gauss)");

   return 0;
}

double che_poww (double A, int b, int *error)
{
	int i = 0;
   int B = 0;
   double vystup=1.0;

   *error = 0;
   B = b;
   //
   if (b<0) B = -b;
   if (b==0) return 1.0;
   if (A==0.0)
   {
	   if (b>0)
      {
		   return 0.0;
      }
	   else
	   {
		   printf ("chyba pri mocneni %f**%d !\n", A, b);
         *error = 1;
//		   exit(111);
		   return 1.0e99;						// cislo vycucane z prstu
	   }
   }
   for (i=0; i<B; i++)
   {
   	vystup *= A;
   }
   if (b<0)
   {
	   vystup=1.0/vystup;
   }

   return vystup;
}

double che_poww_ld (double A, double b, int *error)
{
   int B = 0;

   if (floor(b)==b)
   {
      B=(int) b;
      return che_poww(A,B,error);
   }
   else
   {
   	error = 0;
      if (A<0.0)
      {
		   printf ("chyba pri mocneni %f**%f !\n", A, b);
		 *error = -1;
//		   exit(111);
         return 0.0;						// cislo vycucane z prstu
      }
      else if (A==0.0)
      {
         if (b<0.0)
         {
            printf ("chyba pri mocneni %f**%f !\n", A, b);
	         *error = 1;
//          exit(111);
            return 1.0e99;						// cislo vycucane z prstu
         }
         else
         {
            return 0.0;
         }
      }
      else
      {
         return exp(log(A)*b);
      }
   }
}
/*
void che_cti_param_soubor (char *soubor)
{
   char odpad[1000];
   int i,j;
   FILE *fr;
   TS_che pom_che;

   if ((fr = fopen(soubor, "r")) == NULL)
   {
	   printf ("Nelze otevrit vstupni soubor chemie!\n");
	   exit(0);
   }

   fscanf (fr,"%ld%[^\n]s", &G_prm.pocet_latek, &odpad);//printf ("Kolik latek?");
   if (G_prm.pocet_latek>MAX_POC_LATEK)
   {
	   printf ("Celkovy pocet latek muze byt maximalne %d!\n", MAX_POC_LATEK);
	   exit(121);
   }
   fscanf (fr,"%ld%[^\n]s", &G_prm.pocet_latekvefazi, &odpad);//printf ("Kolik ve vysetrovane fazi?");
   if (G_prm.pocet_latekvefazi>G_prm.pocet_latek)
   {
	   printf ("Pocet latek v jedne fazi je vetsi nez celkovy pocet latek!\n");
	   exit(122);
   }
   fscanf (fr,"%ld%[^\n]s", &G_prm.celkovy_pocet_reakci, &odpad);//printf ("Kolik reakci?");
   fscanf (fr,"%Lf%[^\n]s", &G_prm.T, &odpad);//printf ("Jaka teplota?");
   P_lat = (TS_lat *)malloc( G_prm.pocet_latek*sizeof( TS_lat ) );
   if ( P_lat == NULL )
   {
	   printf ("Malo pameti!\n");
	   exit(0);
   }
   P_che = (TS_che *)malloc( (G_prm.celkovy_pocet_reakci)*sizeof( TS_che ) );
   if ( P_che == NULL )
   {
	   printf ("Malo pameti!\n");
	   exit(0);
   }
   for (i=0; i<G_prm.pocet_latekvefazi; i++)// printf ("%ld", G_prm.pocet_latek);
   {
      fscanf (fr,"%s%[^\n]s", &P_lat[i].nazev,&odpad);//printf ("Nazev latky %d?",i);
	   fscanf (fr,"%ld%[^\n]s", &P_lat[i].Q, &odpad);//printf ("naboj latky %d?",i);
	   fscanf (fr,"%Lf%[^\n]s", &P_lat[i].dGf, &odpad);//printf ("delta G_f latky %d?",i);
	   fscanf (fr,"%Lf%[^\n]s", &P_lat[i].M, &odpad);//printf ("mol. hmotnost M latky %d?",i);
	   fscanf (fr,"%Lf%[^\n]s", &P_lat[i].m0, &odpad);//printf ("poc. molalita latky %d?",i);
   }
   if (G_prm.pocet_latekvefazi<G_prm.pocet_latek)
	   for (i=G_prm.pocet_latekvefazi; i<G_prm.pocet_latek; i++)
	   {
		   fscanf (fr,"%s%[^\n]s", &P_lat[i].nazev,&odpad);//printf ("Nazev latky %d?",i);
		   fscanf (fr,"%Lf%[^\n]s", &P_lat[i].dGf, &odpad);//printf ("delta G_f latky %d?",i);
		   fscanf (fr,"%Lf%[^\n]s", &P_lat[i].M, &odpad);//printf ("mol. hmotnost M latky %d?",i);
		   fscanf (fr,"%Lf%[^\n]s", &P_lat[i].aktivita, &odpad);//printf ("aktivita latky %d?",i);
	   }
   fscanf (fr,"%Lf%[^\n]s", &G_prm.Afi, &odpad);//printf ("Afi?");
   fscanf (fr,"%Lf%[^\n]s", &G_prm.b, &odpad);//printf ("b?");
	G_prm.pocet_reakci_pro_matici = 0;
   for (i=0; i<G_prm.celkovy_pocet_reakci; i++)
   {
		G_prm.pocet_reakci_pro_matici++;
	   for (j=0; j<G_prm.pocet_latek; j++)
		{
		   fscanf (fr,"%ld%[^\n]s", &P_che[i].stech_koef_p[j], &odpad);//printf ("stech. koef. latky %d v reakci %d?",j,i);
      }
	   fscanf (fr,"%ld%[^\n]s", &P_che[i].typ_reakce, &odpad);//printf ("typ reakce c. %d?",i);
	   if (P_che[i].typ_reakce==2)
	   {
		   fscanf (fr,"%Lf%[^\n]s", &P_che[i].K, &odpad);//printf ("bilancni konstanta pro rovnici %d?", i);
		   for (j=0; j<G_prm.pocet_latek; j++)
			{
//    NACITA SE JICH MOC!   UZITECNE JSOU JEN PRO LATKY VE FAZI, OSTATNI SE NEMOHOU UZIT
			   fscanf (fr,"%Lf%[^\n]s", &P_che[i].exponent[j], &odpad);//printf ("vaha latky %d pro bilancni rovnici %d?",j,i);
         }
	   }
	   else if ((P_che[i].typ_reakce==1)||(P_che[i].typ_reakce==3))
	   {
		   fscanf (fr,"%Lf%[^\n]s", &P_che[i].K, &odpad);//printf ("kineticka konstanta pro reakci %d?", i);
		   for (j=0; j<G_prm.pocet_latek; j++)
         {
			   fscanf (fr,"%Lf%[^\n]s", &P_che[i].exponent[j], &odpad);//printf ("mocnina latky %d pro kinetiku reakce %d?",j,i);
         }
         if (P_che[i].typ_reakce==3)
         {
         	G_prm.pocet_reakci_pro_matici--;
         }
	   }
	   else if (P_che[i].typ_reakce==0)
      {
		   fscanf (fr,"%Lf%[^\n]s", &P_che[i].K, &odpad);//printf ("rovnovazna konstanta pro reakci %d?", i);
      }
	   else
	   {
		   printf ("typ reakce muze byt pouze:\n 0 == rovnovazna reakce,\n 1 == kineticka reakce,\n 2 == bilancni rovnice,\n 3 == pomala kineticka reakce pocitana pred rovnovahou");
		   exit (222);
	   }
	   fscanf (fr,"%Lf%[^\n]s", &P_che[i].zeta0, &odpad);//printf ("zeta0 pro reakci %d?", i);
   }
	if (G_prm.celkovy_pocet_reakci>1)
   {
      for (j=0; j<G_prm.celkovy_pocet_reakci; j++)
      {
         for (i=0; i<G_prm.celkovy_pocet_reakci-1; i++)
         {
            if ((P_che[i].typ_reakce == 3) && (P_che[i+1].typ_reakce != 3))
            {
            	pom_che = P_che[i];
               P_che[i] = P_che[i+1];
               P_che[i+1] = pom_che;
            }
         }
      }
	}
   fscanf (fr,"%Lf%[^\n]s", &G_prm.epsilon, &odpad);//printf ("epsilon?");
   fscanf (fr,"%Lf%[^\n]s", &G_prm.omega, &odpad);//printf ("omega?");
   fscanf (fr,"%Lf%[^\n]s", &G_prm.deltaT, &odpad);//printf ("delta t?");
   fscanf (fr,"%ld%[^\n]s", &G_prm.cas_kroku, &odpad);//printf ("pocet cas. kroku?");
   fscanf (fr,"%ld%[^\n]s", &G_prm.vypisy, &odpad);//printf ("uroven vypisu?");

   fclose(fr);
}
*/
double che_m_ (int latka, double *zeta)
{
   int i = 0;
   double mm = 0.0;

   mm=P_lat[latka].m0;
   for (i=0; i<G_prm.pocet_reakci_pro_matici; i++)
   {
	   //printf("\nzeta[%ld]: %lf P_che[%ld].stech_koef_p[%ld]: %lf\n",i,zeta[i],i,latka,P_che[i].stech_koef_p[latka]);
	   mm += zeta[i]*P_che[i].stech_koef_p[latka];
   }
   return mm;
}

double che_m_x (int latka, double *zeta, double krat)
{
   int i = 0;
   double mm = 0.0;

   mm=P_lat[latka].m0;
   for (i=0; i<G_prm.pocet_reakci_pro_matici; i++)
   {
	   mm += krat*zeta[i]*P_che[i].stech_koef_p[latka];
   }
   return mm;
}

double che_m (int latka, double *zeta)
{
   double mm = 0.0;

   mm = che_m_(latka, zeta);
   if (mm<0.0)
   {
//	   printf ("che_m(%ld,zeta)=%Le!\n", latka, mm);
//	   printf ("che_m_x(%ld,zeta,0.1)=%Le!\n", latka, che_m_x(latka, zeta, 0.1));
//	   printf ("che_m_x(%ld,zeta,0.01)=%Le!\n", latka, che_m_x(latka, zeta, 0.01));
//	   printf ("che_m_x(%ld,zeta,0.0)=%Le!\n", latka, che_m_x(latka, zeta, 0.0));
//   	mm = 0.0;
   }
   return mm;
}

double che_I (double *zeta)
{
   double vystup = 0.0;
   int i = 0;

   vystup = 0.0;
   for (i=0; i<G_prm.pocet_latekvefazi; i++)
   {
      vystup += che_m(i, zeta)*P_lat[i].Q*P_lat[i].Q;
   }
	vystup /= 2.0;
   return vystup;
}

double che_dI (int smer)
{
   double vystup = 0.0;
   int i = 0;

   vystup = 0.0;
   for (i=0; i<G_prm.pocet_latekvefazi; i++)
   {
      vystup += P_che[smer].stech_koef_p[i]*P_lat[i].Q*P_lat[i].Q;
   }
   vystup /= 2.0;
   return vystup;
}

double che_gama_ (int i, double *zeta, int *error)
{
   double vystup = 0.0;
   double sqrtlI = 0.0;

   *error = 0;
	if (G_prm.vypisy>4) printf("\nche_gama_:");//xprintf(Msg,"\nche_gama_:");
   if (che_I(zeta)<0.0)
   {
	   printf ("che_I(zeta)=%f!\n", che_I(zeta));
//	   exit (123);
		*error = 1;
   }
   sqrtlI = sqrt(fabs(che_I(zeta)));
   vystup = -P_lat[i].Q*P_lat[i].Q*G_prm.Afi*(sqrtlI/(1.0+G_prm.b*sqrtlI)+2.0/G_prm.b*log(1.0+G_prm.b*sqrtlI));
   vystup = exp(vystup);

	if (G_prm.vypisy>4) xprintf(Msg,"o.k.(che_gama_)");
   return vystup;
}

double che_dgama_ (int i, double *zeta, int smer, int *error)
{
   double vystup = 0.0;
   double sqrtlI = 0.0;
	double pom = 0.0;

   error = 0;
	if (G_prm.vypisy>4) xprintf(Msg,"\nche_dgama_:");
   if (che_dI(smer)==0.0)
   {
	   return 0.0;
   }
   if (che_I(zeta)<0.0)
   {
	   printf ("che_I(zeta)=%f!\n", che_I(zeta));
//	   exit (123);
		*error = 1;
   }
   sqrtlI = sqrt(fabs(che_I(zeta)));
   if (sqrtlI==0.0)
   {
	   xprintf(Msg,"sqrtlI = 0.0, posouvam ve smeru: ");
		pom = zeta[smer];
      zeta[smer]+=DROBNY_POSUN;						// cislo vycucane z prstu
	   sqrtlI = sqrt(fabs(che_I(zeta)));
      zeta[smer]=pom;
	   xprintf(Msg,"sqrtlI = %f\n", sqrtlI);
      if (che_I(zeta)<0.0)
      {
         printf ("che_I(zeta)=%f, posouvam proti smeru: ", che_I(zeta));
         zeta[smer]-=DROBNY_POSUN;						// cislo vycucane z prstu
         sqrtlI = sqrt(fabs(che_I(zeta)));
         zeta[smer]=pom;
         xprintf(Msg,"sqrtlI = %f\n", sqrtlI);
      }
   }
   vystup = -P_lat[i].Q*P_lat[i].Q*G_prm.Afi*(sqrtlI/(1.0+G_prm.b*sqrtlI)+2.0/G_prm.b*log(1.0+G_prm.b*sqrtlI));
   vystup = exp(vystup);
   vystup*= -P_lat[i].Q*P_lat[i].Q*G_prm.Afi;
//   vystup*= (-G_prm.b*sqrtlI/(1.0+G_prm.b*sqrtlI)+3.0)/(1.0+G_prm.b*sqrtlI);
   vystup*= (3.0+2.0*G_prm.b*sqrtlI)/(1.0+G_prm.b*sqrtlI)/(1.0+G_prm.b*sqrtlI);
   if (sqrtlI==0.0)
   {
	   xprintf(Msg,"sqrtlI = 0.0: ");
//	   exit(133);
		vystup *= 1.0e56;						// cislo vycucane z prstu
	   xprintf(Msg,"vystup = %f\n", vystup);
   }
   else
   {
   	vystup/= 2.0*sqrtlI;
   }
   vystup*= che_dI(smer);
	if (G_prm.vypisy>4) xprintf(Msg,"o.k.(che_dgama_)");

   return vystup;
}

double che_K1_( double *zeta, int rce, int vnoreni)
{
   int i = 0;
   double k1 = 0.0;
   int chyba = 0;
   double pom = 0.0;

   if (vnoreni>2*PUL_POCTU_VNORENI)									// cislo vycucane z prstu
   {
      xprintf(Msg,"k1 se moc spatne pocita, koncim!");
   	exit (222);
   }
	if (G_prm.vypisy>4) xprintf(Msg,"\n (che_K1_)");
   if (P_che[rce].typ_reakce==2)
   {
//		if (G_prm.vypisy>4) printf("-reakce 2");
	   k1=0.0;
	   for (i=0; i<G_prm.pocet_latekvefazi; i++)
      {
		   k1+=che_m(i,zeta)*P_che[rce].exponent[i];
      }
	   return k1;
   }
   else if (P_che[rce].typ_reakce==1)
   {
//		if (G_prm.vypisy>4) printf("-reakce 1");
	   k1=1.0;
	   for (i=0; i<G_prm.pocet_latekvefazi; i++)
      {
		   k1*=che_poww_ld(che_m(i,zeta),P_che[rce].exponent[i],&chyba);
         if (chyba > 0)
         {
         	// nula na zaporne cislo nebo zaporne cislo na necele cislo!
            pom = zeta[rce];
            if (vnoreni < PUL_POCTU_VNORENI) 		// cislo vycucane z prstu
            {
               xprintf(Msg,"k1 se spatne pocita, posouvam ve smeru: ");
               zeta[rce]+=DROBNY_POSUN;						// cislo vycucane z prstu
            }
            else if (vnoreni == PUL_POCTU_VNORENI)	// cislo vycucane z prstu
            {
               xprintf(Msg,"k1 se spatne pocita, posouvam poprve proti smeru: ");
               zeta[rce]-=(PUL_POCTU_VNORENI+1)*DROBNY_POSUN;	// cislo vycucane z prstu
            }
            else
            {
               xprintf(Msg,"k1 se spatne pocita, posouvam proti smeru: ");
               zeta[rce]-=DROBNY_POSUN;						// cislo vycucane z prstu
            }
            k1 = che_K1_(zeta, rce, vnoreni+1);
            zeta[rce]=pom;
            return k1;
         }
      }
//		if (G_prm.vypisy>4) printf("<%Le/%Le>", zeta[rce],k1);
      if (k1 == 0.0)
      {
      	if (zeta[rce]==0.0)
         {
         	k1 = 0.0;
         }
         else
         {
         	k1 = zeta[rce]*1e99;     // cislo vycucane z prstu
         }
      }
      else
      {
		   k1 = zeta[rce]/k1;
      }
	   return k1;
   }
   else
   {
//		if (G_prm.vypisy>4) printf("-jina reakce");
	   k1=1.0;
	   for (i=0; i<G_prm.pocet_latekvefazi; i++)
      {
		   k1*=che_poww(che_m(i,zeta),P_che[rce].stech_koef_p[i],&chyba);
         if (chyba > 0)
         {
         	// nula na zaporne cislo!
            pom = zeta[rce];
            if (vnoreni < PUL_POCTU_VNORENI)			// cislo vycucane z prstu
            {
               xprintf(Msg,"k1 se spatne pocita, posouvam ve smeru: ");
               zeta[rce]+=DROBNY_POSUN;						// cislo vycucane z prstu
            }
            else if (vnoreni == PUL_POCTU_VNORENI)	// cislo vycucane z prstu
            {
               xprintf(Msg,"k1 se spatne pocita, posouvam poprve proti smeru: ");
               zeta[rce]-=(PUL_POCTU_VNORENI+1)*DROBNY_POSUN;// cislo vycucane z prstu
            }
            else
            {
               xprintf(Msg,"k1 se spatne pocita, posouvam proti smeru: ");
               zeta[rce]-=DROBNY_POSUN;						// cislo vycucane z prstu
            }
            k1 = che_K1_(zeta, rce, vnoreni+1);
            zeta[rce]=pom;
            return k1;
         }
      }
	   return k1;
   }
}

double che_K2_( double *zeta, int rce, int vnoreni)
{
   int i = 0;
   double k2 = 0.0;
   int chyba = 0;
   int chybicka = 0;
   double pom = 0.0;

   if (vnoreni>2*PUL_POCTU_VNORENI)					// cislo vycucane z prstu
   {
      xprintf(Msg,"k2 se moc spatne pocita, koncim!");
   	exit (222);
   }
	if (G_prm.vypisy>4) xprintf(Msg,"\n (che_K2_)");
   if (P_che[rce].typ_reakce!=0)
   {
   	return 1.0;
   }
   k2=1.0;
   for (i=0; i<G_prm.pocet_latekvefazi; i++)
   {
		k2*=che_poww(che_gama_(i,zeta,&chybicka),P_che[rce].stech_koef_p[i],&chyba);
      if (chybicka > 0)
		{
      	// zaporna iontova sila!
         xprintf(Msg,"Zaporna iontova sila - ");
         chyba += 1;
      }
      if (chyba > 0)
      {
         // nula na zaporne cislo!
         pom = zeta[rce];
         if (vnoreni < PUL_POCTU_VNORENI)				// cislo vycucane z prstu
         {
            xprintf(Msg,"k2 se spatne pocita, posouvam ve smeru: ");
            zeta[rce]+=DROBNY_POSUN;						// cislo vycucane z prstu
         }
         else if (vnoreni == PUL_POCTU_VNORENI)	// cislo vycucane z prstu
         {
            xprintf(Msg,"k2 se spatne pocita, posouvam poprve proti smeru: ");
            zeta[rce]-=(PUL_POCTU_VNORENI+1)*DROBNY_POSUN;		// cislo vycucane z prstu
         }
         else
         {
            xprintf(Msg,"k2 se spatne pocita, posouvam proti smeru: ");
            zeta[rce]-=DROBNY_POSUN;						// cislo vycucane z prstu
         }
         k2 = che_K2_(zeta, rce, vnoreni+1);
         zeta[rce]=pom;
         return k2;
      }
   }
   return k2;
}

double che_lnKT0(int rce)
{
   int i = 0;
   double kk = 0.0;

   kk=0.0;
   for (i=0; i<G_prm.pocet_latek; i++)
   {
      kk+=P_che[rce].stech_koef_p[i]*P_lat[i].dGf;
   }
   kk/=-1.0*_R*G_prm.TGf;

   return kk;
}

double che_dH(int rce)
{
   int i = 0;
   double hh = 0.0;

   hh=0.0;
   for (i=0; i<G_prm.pocet_latek; i++)
   {
      hh+=P_che[rce].stech_koef_p[i]*P_lat[i].dHf;
   }

   return hh;
}

double che_K_(int rce)
{
   int i = 0;
   double kk = 0.0;
   int chyba = 0;

	if (G_prm.vypisy>4) xprintf(Msg,"\n (che_K_)");
   if (P_che[rce].typ_reakce==2)
   {
	   kk=P_che[rce].K;
	   return kk;
   }
   else if (P_che[rce].typ_reakce==1)
   {
      kk=1.0;
      for (i=0; i<G_prm.pocet_latekvefazi; i++)
      {
		 kk *= che_poww_ld(P_lat[i].m0,fabs(P_che[rce].exponent[i]),&chyba);
//printf("kk*=%Le",che_poww_ld(P_lat[i].m0,fabsl(P_che[rce].exponent[i]),chyba));
         if (chyba > 0)
         {
            // nula na zaporne cislo nebo zaporne cislo na necely exponent!
			xprintf(Msg,"K se moc spatne pocita, koncim!");
            exit (222);
         }
      }
      if (kk > ACCURACY)
      {
	      kk=P_che[rce].K*G_prm.deltaT;
      }
//printf("kk=%Le",kk);
	   return kk;
   }
   else
   {
	   if (P_che[rce].K>0.0)
  		   kk = P_che[rce].K;
	   else
  	   {
		   kk=exp(che_lnKT0(rce)-che_dH(rce)*_R*(1.0/G_prm.T-1.0/G_prm.TGf));
	   }

	   if (G_prm.pocet_latekvefazi<G_prm.pocet_latek)
      {
		   for (i=G_prm.pocet_latekvefazi; i<G_prm.pocet_latek; i++)
         {
			   kk/=che_poww(P_lat[i].aktivita,P_che[rce].stech_koef_p[i],&chyba);
				if (chyba > 0)
            {
               // nula na zaporne cislo!
               xprintf(Msg,"K se moc spatne pocita, koncim!");
               exit (222);
            }
         }
      }
	   return kk;
   }
}

double che_dK1_(double *zeta, int rce, int smer, int vnoreni)
{
   int i = 0;
   int j = 0;
   double dk1 = 0.0;
   double pomm = 0.0;
   int chyba = 0;
   double pom = 0.0;

   if (vnoreni>2*PUL_POCTU_VNORENI)					// cislo vycucane z prstu
   {
      xprintf(Msg,"dk1 se moc spatne pocita, koncim!");
   	exit (222);
   }
	if (G_prm.vypisy>4) xprintf(Msg,"\n (che_dK1_)");
   if (P_che[rce].typ_reakce==2)
   {
   	dk1 = 0.0;
		for (i=0; i<G_prm.pocet_latekvefazi; i++)
		{
      	dk1+=P_che[rce].exponent[i]*P_che[smer].stech_koef_p[i];
      }
		return dk1;
   }
   else if (P_che[rce].typ_reakce==1)
   {
	   dk1 = 0.0;
	   if ( rce == smer )
	   {
		   dk1=1.0;
		   for (i=0; i<G_prm.pocet_latekvefazi; i++)
         {
			  dk1*=che_poww_ld(che_m(i,zeta),-1.0*P_che[rce].exponent[i],&chyba);
            if (chyba > 0)
            {
               // nula na zaporne cislo nebo zaporne cislo na necele cislo!
               pom = zeta[rce];
               if (vnoreni < PUL_POCTU_VNORENI)			// cislo vycucane z prstu
               {
                  xprintf(Msg,"dk1 se spatne pocita, posouvam ve smeru: ");
                  zeta[rce]+=DROBNY_POSUN;						// cislo vycucane z prstu
               }
               else if (vnoreni == PUL_POCTU_VNORENI)	// cislo vycucane z prstu
               {
                  xprintf(Msg,"dk1 se spatne pocita, posouvam poprve proti smeru: ");
                  zeta[rce]-=(PUL_POCTU_VNORENI+1)*DROBNY_POSUN;	// cislo vycucane z prstu
               }
               else
               {
                  xprintf(Msg,"dk1 se spatne pocita, posouvam proti smeru: ");
                  zeta[rce]-=DROBNY_POSUN;						// cislo vycucane z prstu
               }
               dk1 = che_dK1_(zeta, rce, smer, vnoreni+1);
               zeta[rce]=pom;
               return dk1;
            }
         }
	   }
	   for (i=0; i<G_prm.pocet_latekvefazi; i++)
	   {
		   pomm = 1.0;
		   for (j=0; j<G_prm.pocet_latekvefazi; j++)
         {
			   if (j!=i)
            {
				   pomm*=che_poww_ld(che_m(j,zeta),-1.0*P_che[rce].exponent[j],&chyba);
               if (chyba > 0)
               {
                  // nula na zaporne cislo nebo zaporne cislo na necele cislo!
                  pom = zeta[rce];
                  if (vnoreni < PUL_POCTU_VNORENI)			// cislo vycucane z prstu
                  {
                     xprintf(Msg,"dk1 se spatne pocita, posouvam ve smeru: ");
                     zeta[rce]+=DROBNY_POSUN;						// cislo vycucane z prstu
                  }
                  else if (vnoreni == PUL_POCTU_VNORENI)	// cislo vycucane z prstu
                  {
                     xprintf(Msg,"dk1 se spatne pocita, posouvam poprve proti smeru: ");
                     zeta[rce]-=(PUL_POCTU_VNORENI+1)*DROBNY_POSUN;	// cislo vycucane z prstu
                  }
                  else
                  {
                     xprintf(Msg,"dk1 se spatne pocita, posouvam proti smeru: ");
                     zeta[rce]-=DROBNY_POSUN;						// cislo vycucane z prstu
                  }
                  dk1 = che_dK1_(zeta, rce, smer, vnoreni+1);
                  zeta[rce]=pom;
                  return dk1;
               }
            }
         }
		   if (P_che[rce].stech_koef_p[i]!=0)
         {
			  //double docasna = -1.0*P_che[rce].exponent[j]-1.0;  if(docasna == 0.0) docasna = 0.0; 
			  pomm*=che_poww_ld(che_m(i,zeta),-1.0*P_che[rce].exponent[i]-1.0,&chyba)*(-1.0*P_che[rce].exponent[i])*(P_che[smer].stech_koef_p[i]);
            if (chyba > 0)
            {
               // nula na zaporne cislo nebo zaporne cislo na necele cislo!
               pom = zeta[rce];
               if (vnoreni < PUL_POCTU_VNORENI)			// cislo vycucane z prstu
               {
                  xprintf(Msg,"dk1 se spatne pocita, posouvam ve smeru: ");
                  zeta[rce]+=DROBNY_POSUN;						// cislo vycucane z prstu
               }
               else if (vnoreni == PUL_POCTU_VNORENI)	// cislo vycucane z prstu
               {
                  xprintf(Msg,"dk1 se spatne pocita, posouvam poprve proti smeru: ");
                  zeta[rce]-=(PUL_POCTU_VNORENI+1)*DROBNY_POSUN;	// cislo vycucane z prstu
               }
               else
               {
                  xprintf(Msg,"dk1 se spatne pocita, posouvam proti smeru: ");
                  zeta[rce]-=DROBNY_POSUN;						// cislo vycucane z prstu
               }
               dk1 = che_dK1_(zeta, rce, smer, vnoreni+1);
               zeta[rce]=pom;
               return dk1;
            }
         }
		   else
         {
			   pomm = 0.0;
         }
		   dk1 += zeta[rce]*pomm;
	   }
	   return dk1;
   }
   else // typ_reakce == 0
   {
	   dk1=0.0;
	   for (i=0; i<G_prm.pocet_latekvefazi; i++)
	   {
		   pomm = 1.0;
		   for (j=0; j<G_prm.pocet_latekvefazi; j++)
         {
			   if (j!=i)
            {
//printf("%Le ** %ld\n",che_m(j,zeta),P_che[rce].stech_koef_p[j]);
				   pomm*=che_poww(che_m(j,zeta),P_che[rce].stech_koef_p[j],&chyba);
               if (chyba > 0)
               {
                  // nula na zaporne cislo!
                  pom = zeta[rce];
                  if (vnoreni < PUL_POCTU_VNORENI)		// cislo vycucane z prstu
                  {
                     xprintf(Msg,"dk1 se spatne pocita, posouvam ve smeru: ");
                     zeta[rce]+=DROBNY_POSUN;						// cislo vycucane z prstu
                  }
                  else if (vnoreni == PUL_POCTU_VNORENI)	// cislo vycucane z prstu
                  {
                     xprintf(Msg,"dk1 se spatne pocita, posouvam poprve proti smeru: ");
                     zeta[rce]-=(PUL_POCTU_VNORENI+1)*DROBNY_POSUN;	// cislo vycucane z prstu
                  }
                  else
                  {
                     xprintf(Msg,"dk1 se spatne pocita, posouvam proti smeru: ");
                     zeta[rce]-=DROBNY_POSUN;						// cislo vycucane z prstu
                  }
                  dk1 = che_dK1_(zeta, rce, smer, vnoreni+1);
                  zeta[rce]=pom;
                  return dk1;
               }
            }
         }
		   if (P_che[rce].stech_koef_p[i]!=0)
         {
			   pomm*=che_poww(che_m(i,zeta),P_che[rce].stech_koef_p[i]-1,&chyba)*(P_che[rce].stech_koef_p[i])*(P_che[smer].stech_koef_p[i]);
            if (chyba > 0)
            {
               // nula na zaporne cislo!
               pom = zeta[rce];
               if (vnoreni < PUL_POCTU_VNORENI)			// cislo vycucane z prstu
               {
                  xprintf(Msg,"dk1 se spatne pocita, posouvam ve smeru: ");
                  zeta[rce]+=DROBNY_POSUN;						// cislo vycucane z prstu
               }
               else if (vnoreni == PUL_POCTU_VNORENI)						// cislo vycucane z prstu
               {
                  xprintf(Msg,"dk1 se spatne pocita, posouvam poprve proti smeru: ");
                  zeta[rce]-=(PUL_POCTU_VNORENI+1)*DROBNY_POSUN;						// cislo vycucane z prstu
               }
               else
               {
                  xprintf(Msg,"dk1 se spatne pocita, posouvam proti smeru: ");
                  zeta[rce]-=DROBNY_POSUN;						// cislo vycucane z prstu
               }
               dk1 = che_dK1_(zeta, rce, smer, vnoreni+1);
               zeta[rce]=pom;
               return dk1;
            }
         }
		   else
         {
			   pomm = 0.0;
         }
		   dk1+=pomm;
	   }
	   return dk1;
   }
}

double che_dK2_(double *zeta, int rce, int smer, int vnoreni)
{
   int i = 0;
   int j = 0;
   double dk2 = 0.0;
   double pomm = 0.0;
   int chyba = 0;
   int chybicka = 0;
   int chybicka2 = 0;
   double pom = 0.0;

   if (vnoreni>2*PUL_POCTU_VNORENI)									// cislo vycucane z prstu
   {
      xprintf(Msg,"dk2 se moc spatne pocita, koncim!");
   	exit (222);
   }
	if (G_prm.vypisy>4) xprintf(Msg,"\n (che_dK2_)");
   if (P_che[rce].typ_reakce) return 0.0;
   dk2=0.0;
   for (i=0; i<G_prm.pocet_latekvefazi; i++)
   {
	   pomm = 1.0;
	   for (j=0; j<G_prm.pocet_latekvefazi; j++)
      {
		   if (j!=i)
         {
			   pomm*=che_poww(che_gama_(j,zeta,&chybicka),P_che[rce].stech_koef_p[j],&chyba);
            if (chybicka > 0)
            {
               // zaporna iontova sila!
               xprintf(Msg,"Zaporna iontova sila - ");
               chyba += 1;
            }
            if (chyba > 0)
            {
               // nula na zaporne cislo!
               pom = zeta[rce];
               if (vnoreni < PUL_POCTU_VNORENI)								// cislo vycucane z prstu
               {
                  xprintf(Msg,"dk2 se spatne pocita, posouvam ve smeru: ");
                  zeta[rce]+=DROBNY_POSUN;						// cislo vycucane z prstu
               }
               else if (vnoreni == PUL_POCTU_VNORENI)						// cislo vycucane z prstu
               {
                  xprintf(Msg,"dk2 se spatne pocita, posouvam poprve proti smeru: ");
                  zeta[rce]-=(PUL_POCTU_VNORENI+1)*DROBNY_POSUN;						// cislo vycucane z prstu
               }
               else
               {
                  xprintf(Msg,"dk2 se spatne pocita, posouvam proti smeru: ");
                  zeta[rce]-=DROBNY_POSUN;						// cislo vycucane z prstu
               }
               dk2 = che_dK2_(zeta, rce, smer, vnoreni+1);
               zeta[rce]=pom;
               return dk2;
            }
         }
      }
      if (P_che[rce].stech_koef_p[i]!=0)
      {
	      dk2+=pomm*P_che[rce].stech_koef_p[i]*che_poww(che_gama_(i,zeta,&chybicka),P_che[rce].stech_koef_p[i]-1,&chyba)*che_dgama_(i,zeta,smer,&chybicka2);
         if ((chybicka + chybicka2) > 0)
         {
            // zaporna iontova sila!
            xprintf(Msg,"Zaporna iontova sila - ");
            chyba += 1;
         }
         if (chyba > 0)
         {
         	// nula na zaporne cislo!
            pom = zeta[rce];
            if (vnoreni < PUL_POCTU_VNORENI)								// cislo vycucane z prstu
            {
               xprintf(Msg,"dk2 se spatne pocita, posouvam ve smeru: ");
               zeta[rce]+=DROBNY_POSUN;						// cislo vycucane z prstu
            }
            else if (vnoreni == PUL_POCTU_VNORENI)						// cislo vycucane z prstu
            {
               xprintf(Msg,"dk2 se spatne pocita, posouvam poprve proti smeru: ");
               zeta[rce]-=(PUL_POCTU_VNORENI+1)*DROBNY_POSUN;						// cislo vycucane z prstu
            }
            else
            {
               xprintf(Msg,"dk2 se spatne pocita, posouvam proti smeru: ");
               zeta[rce]-=DROBNY_POSUN;						// cislo vycucane z prstu
            }
            dk2 = che_dK2_(zeta, rce, smer, vnoreni+1);
            zeta[rce]=pom;
            return dk2;
         }
      }
   }
   return dk2;
}

double che_hodnota_(double *zeta, int rce)
{
	double vystup = 0.0;

	if (G_prm.vypisy>4) xprintf(Msg,"\nche_hodnota_:");
	vystup = che_K1_(zeta,rce,0)*che_K2_(zeta,rce,0);
	if (G_prm.vypisy>4) xprintf(Msg,"o.k. (che_hodnota_)");

   return vystup;
}

double che_derivace_(double *zeta, int rce, int smer)
{
	double vystup = 0.0;

	if (G_prm.vypisy>4) xprintf(Msg,"\nche_derivace_:");
   vystup = che_K1_(zeta,rce,0)*che_dK2_(zeta,rce,smer,0)+che_K2_(zeta,rce,0)*che_dK1_(zeta,rce,smer,0);
	if (G_prm.vypisy>4) xprintf(Msg,"o.k. (che_derivace_)");

   return vystup;
}

void che_hodnoty(double *zeta, double *hodnota, double *skala)
{
   int i = 0;

   for (i=0;i<G_prm.pocet_reakci_pro_matici;i++)
   {
	   hodnota[i]=che_hodnota_(zeta,i)*skala[i];
   }
}

void che_Jakobi(double *zeta, double *J, double *skala)
{
   int i = 0;
   int j =0;

	if (G_prm.vypisy>4) xprintf(Msg,"\nche_Jakobi: ");
   for (i=0;i<G_prm.pocet_reakci_pro_matici;i++)
   {
	   for (j=0;j<G_prm.pocet_reakci_pro_matici;j++)
      {
		   J[i*G_prm.pocet_reakci_pro_matici+j]=che_derivace_(zeta,i,j)*skala[i];
      }
   }
	if (G_prm.vypisy>4) xprintf(Msg,"o.k. (Jakobi)");
}

double che_abs_norma(double *x)
{
   int i = 0;
   double vysl = 0.0;

   for (i=0;i<G_prm.pocet_reakci_pro_matici;i++)
   {
	   vysl+=x[i]*x[i];
   }

   return sqrt(vysl);
}

double che_rel_norma(double *x, double *K)
{
   int i = 0;

   double vysl = 0.0;
   for (i=0;i<G_prm.pocet_reakci_pro_matici;i++)
	{
      vysl+=x[i]*x[i]/K[i]/K[i];
   }

   return sqrt(vysl);
}

double che_norma (double *x, double *K)
{
	double norma = 0.0;
   FILE *fw;

	if (G_prm.abs_norma) norma=che_abs_norma(x);
   else norma=che_rel_norma(x,K);
   if (G_prm.vypisy>3)
   {
   	fw = fopen(vystupni_soubor, "a");
   	fprintf (fw, "\n>  norma = %f", norma);
      fclose (fw);
   }
   return norma;
}

void che_odecti (double *x, double *y, double *z)
{
   int i = 0;

	if (G_prm.vypisy>4) xprintf(Msg,"\nche_odecti: ");
   for (i=0;i<G_prm.pocet_reakci_pro_matici;i++)
	{
      z[i]=x[i]-y[i];
   }
	if (G_prm.vypisy>4) xprintf(Msg,"o.k. (che_odecti)");
}

void che_nasob_ld (double x, double *y, double *z, int delka)
{
   int i = 0;

	if (G_prm.vypisy>4) xprintf(Msg,"\nche_nasob_ld: ");
   for (i=0;i<delka;i++)
	{
      z[i]=x*y[i];
   }
	if (G_prm.vypisy>4) xprintf(Msg,"o.k. (che_nasob_ld)");
}

void che_kombinuj4_ld (double x1, double *y1, double x2, double *y2, double x3, double *y3, double x4, double *y4, double *z, int delka)
{
   int i = 0;

	if (G_prm.vypisy>4) xprintf(Msg,"\nche_kombinuj4_ld: ");
   for (i=0;i<delka;i++)
	{
      z[i]=x1*y1[i]+x2*y2[i]+x3*y3[i]+x4*y4[i];
   }
	if (G_prm.vypisy>4) xprintf(Msg,"o.k. (che_kombinuj4_ld)");
}

void che_nuluj_ld (double *z, int delka)
{
   int i = 0;

	if (G_prm.vypisy>4) xprintf(Msg,"\nche_nuluj_ld: ");
   for (i=0;i<delka;i++)
	{
      z[i]=0.0;
   }
	if (G_prm.vypisy>4) xprintf(Msg,"o.k. (che_nuluj_ld)");
}

void che_kopiruj (double *y, double *z)
{
   int i = 0;

	if (G_prm.vypisy>4) xprintf(Msg,"\nkopiruj:");
   for (i=0;i<G_prm.pocet_reakci_pro_matici;i++)
   {
      z[i]=y[i];
   }
	if (G_prm.vypisy>4) xprintf(Msg,"o.k. (kopiruj)");
}

void che_zgaussproprg ( double *prstrana, int *hprvky, double *vysl )
{
   int ij = 0;

   for ( ij = 0; ij < G_prm.pocet_reakci_pro_matici; ++ij )
   {
	   vysl[ij]=prstrana[hprvky[ij]];
   }
}

void che_napis_soubor_ld(FILE *fw, double *vektor, int delka)
{
   int i;
   fprintf(fw,"(");
   for (i=0;i<delka;i++)
   {
	   fprintf(fw,"%f",vektor[i]);
	   if (i<(delka-1))
      {
		   fprintf(fw,",");
      }
   }
   fprintf(fw,")");
}

void che_napismatici_soubor (char *soubor, double *matice, double *prstr)
{
   int i = 0;
   int j = 0;
   FILE *fw;

   fw = fopen(soubor, "a");
   for ( i=0; i<G_prm.pocet_reakci_pro_matici; ++i )
   {
      fprintf (fw,"\n");
      for ( j=0; j<G_prm.pocet_reakci_pro_matici; ++j )
      {
	      fprintf (fw,"%f ", matice[i*G_prm.pocet_reakci_pro_matici+j]);
      }
      fprintf (fw,"| %f", prstr[i]);
   }
   fclose(fw);
}

int che_odecti_s_korekci_ld(double *x, double *y, double *z, int delka)
{
   int i = 0;
   int j = 0;
   int pruchodu = 0;
   int problem = 0;
   double *z0 = NULL;
   FILE *fw = NULL;

	if (G_prm.vypisy>4) xprintf(Msg,"\nche_odecti_s_korekci_ld:");
   for (i=0;i<G_prm.pocet_latekvefazi;i++)
   {
      if (che_m(i,x) < 0.0)
      {
         printf ("m(%d,...)=%f je zaporne!\n",i,che_m(i,x));
         exit(112);
      }
      if (G_prm.vypisy>3) if (che_m(i,x) == 0.0)
      {
         printf ("m(%d,...)=%f je nulove!\n",i,che_m(i,x));
      }
	}
	z0 = (double *)malloc( delka*sizeof( double ));
   if ( z0 == NULL )
   {
   	printf ("Malo pameti!\n");
   	exit(0);
   }for(j = 0; j < delka; j++){ z0[j] = 0.0; }

   che_odecti(x,y,z0);
   pruchodu = 0;
   do
   {
   	pruchodu++;
      if (pruchodu>MAX_POC_VNITR_CYK)						   // cislo vycucane z prstu
      {
         free (z0);
	 z0 = NULL;
         return 1;
      }
		if (G_prm.vypisy>4)
      {
		   fw = fopen(vystupni_soubor, "a");
      	fprintf(fw,"\n      %d. VNITRNI CYKLUS", pruchodu);
			fprintf(fw,"\nx=");
			che_napis_soubor_ld (fw, x, delka);
			fprintf(fw,"\ny=");
			che_napis_soubor_ld (fw, y, delka);
			fprintf(fw,"\nz0=");
			che_napis_soubor_ld (fw, z0, delka);
         fclose(fw);
   		if (G_prm.vypisy>6)
         {
            fw = fopen(vystupni_soubor, "a");
            fprintf(fw,"\nmolality:");
            for (i=0;i<G_prm.pocet_latekvefazi;i++)
            {
               fprintf(fw,"\t%f",che_m(i,z0));
            }
            fclose (fw);
         }
      }
	   problem = 0;
	   for (i=0;i<G_prm.pocet_latekvefazi;i++)
	   {
	      if (che_m(i,z0) < 0.0)
         {
				if (G_prm.vypisy>4)
            {
					fw = fopen(vystupni_soubor, "a");
            	fprintf(fw, "\nproblem: m(%d)=%f",i,che_m(i,z0));
               fclose (fw);
            }
         	problem = 1;
         }
      }
      if (problem > 0)
      {
      	che_nasob_ld(0.5, y, y, delka);
		   che_odecti(x,y,z0);
      }
   }
   while (problem > 0);
	che_kopiruj (z0, z);
	if (G_prm.vypisy>4) xprintf(Msg,"O.K.(che_odecti_s_korekci_ld)");

   free(z0);
   z0 = NULL;
   return 0;
}

void che_napis_stav_soubor(char *soubor, double *zeta, double *K, double *matice, double *prstr)
{
   FILE *fw = NULL;

   fw = fopen(soubor, "a");
   fprintf (fw, "\n..zeta=");
   che_napis_soubor_ld(fw, zeta, G_prm.pocet_reakci_pro_matici);
   fprintf (fw, "\n.....K=");
   che_napis_soubor_ld(fw, K, G_prm.pocet_reakci_pro_matici);
   fprintf (fw, "\n.....J=");
   fclose(fw);
   che_napismatici_soubor (soubor, matice, prstr);
}

void che_srovnej (double *KK, double *skala)
{
	int i = 0;

   if (G_prm.skaluj_matici > 0)
   {
      for ( i=0; i<G_prm.pocet_reakci_pro_matici; i++)
      {
         if (KK[i] != 0.0)
         {
            skala[i] = 1.0/KK[i];
            KK[i]=1.0;
         }
      }
   }
   else
   {
      for ( i=0; i<G_prm.pocet_reakci_pro_matici; i++)
      {
         skala[i] = 1.0;
      }
   }
}

void che_presun_poc_p_(void)
{
   int i = 0;

   for (i=0; i<G_prm.pocet_latekvefazi; i++)
   {
	   P_lat[i].m0=P_lat[i].m;
	   P_lat[i].m0_sorb=P_lat[i].m_sorb;
   }
}

double che_osklivost(double *zeta0, int *zapornych, int *nulovych, int *nejhorsi)
{
	int i = 0;
   double pom = 0.0;
   double hodnota = 0.0;
   double vysledek = 0.0;

	if (G_prm.vypisy>4) xprintf(Msg,"\nche_osklivost: ");
	vysledek = 0.0;
   *nejhorsi = -1;
   hodnota = 1.0;
   *zapornych = 0;
   *nulovych = 0;
	//printf("\nnulovych %d\n",*nulovych);
   for ( i=0; i<G_prm.pocet_latekvefazi; i++)
   {
	//(*nulovych)++; //JENOM POKUS JESTLI TA INKREMENTACE FUNGUJE
	pom = che_m(i,zeta0);
	//printf("\npom = che_m(i,zeta0), latka %d, pom je %Lf\n",i,pom);
	if (pom <= 0.0)
	  {
		//printf("\npom je mensi nebo rovne nule\n");
		if (pom == 0.0)
		 {
			//printf("\npom je rovne nule\n");
			(*nulovych)++;
		 }
		 else
		 {
			//printf("\npom je mensi nez nula\n");
			(*zapornych)++;
			vysledek -= pom;
		 }
			if (pom<hodnota)
		 {
			//printf("\npom je mensi nez hodnota\n");
			*nejhorsi = i;
			hodnota = pom;
		 }
	  }
	}

	if (vysledek>0.0)
   {
	//printf("\nvysledek je vetsi nez nula\n");
	vysledek += 1.0;
   }
   else
   {
		//printf("\nvysledek je mensi nebo rovny nule\n");
		//neznama promenna nulovych
		vysledek = (1.0 * (*nulovych)) / G_prm.pocet_latekvefazi;
		//vysledek = 0.9;
		//if(G_prm.pocet_latekvefazi != 0) printf("\nlatek ve fazi %d, nulovych %d\n",G_prm.pocet_latekvefazi,*nulovych);
   }
	//printf("\nkonec che_osklivost\n");
	if (G_prm.vypisy>4) xprintf(Msg,"o.k.(che_osklivost = %f)", vysledek);
   return vysledek;
}

int che_urci_zeta0(void)
{
	int i = 0;
   int j = 0;
   int k = 0;
   int x = 0;
   int y = 0;
   double *zeta0 = NULL;
   double *zeta = NULL;
   int zapornych0 = 0;
   int nulovych0 = 0;
   int nejhorsi0 = 0;
   double osklivost0 = 0.0;
   double osklivost = 0.0;
   int chyba = 0;
   int vratit = 0;

   vratit = 0;
	zeta0 = (double *)malloc( G_prm.pocet_reakci_pro_matici*sizeof( double ));
   if ( zeta0 == NULL )
   {
   	printf ("Malo pameti!\n");
   	exit(0);
   }for(j = 0; j < G_prm.pocet_reakci_pro_matici; j++){ zeta0[j] = 0.0; }

	zeta = (double *)malloc( G_prm.pocet_reakci_pro_matici*sizeof( double ));
   if ( zeta == NULL )
   {
   	printf ("Malo pameti!\n");
   	exit(0);
   }for(j = 0; j < G_prm.pocet_reakci_pro_matici; j++){ zeta[j] = 0.0; }
   for ( i=0; i<G_prm.pocet_reakci_pro_matici; i++)
   {
   	zeta0[i]=0.0;
   }
   osklivost0=che_osklivost(zeta0, &zapornych0, &nulovych0, &nejhorsi0);
   if (osklivost0 <= 1.0)
   {
      if (osklivost0 != 0.0)
      {
         for ( i=1; i<che_poww(3.0,G_prm.pocet_reakci_pro_matici,&chyba); i++)
         {
            x = i;
            for ( j=0; j<G_prm.pocet_reakci_pro_matici; j++)
            {
               y=fmod(x,3.0);
               x /=3;
               if ((P_che[j].typ_reakce == 1)&& (y ==2))
               {
                  continue;
               }
               switch (y)
               {
                  case 0: zeta[j] =0.0; break;
                  case 1: zeta[j] =DROBNY_POSUN; break;    //cislo vycucane z prstu
                  case 2: zeta[j] =-DROBNY_POSUN; break;    //cislo vycucane z prstu
               }
            }
			osklivost=che_osklivost (zeta, &zapornych0, &nulovych0, &nejhorsi0);
            if (osklivost < osklivost0)
            {
   //printf("\n%f<%Le ", osklivost, osklivost0);
               osklivost0 = osklivost;
               for ( k=0; k<G_prm.pocet_reakci_pro_matici; k++)
               {
                  zeta0[k]=zeta[k];
               }
            }
            if (osklivost == 0.0) break;
         }
      }
   }
   else
   {
      vratit = 1;
   }
   for ( i=0; i<G_prm.pocet_reakci_pro_matici; i++)
   {
   	P_che[i].zeta0 = zeta0[i];
//   	printf("\nzeta0[%d]=%Le", i,zeta0[i]);
   }
   if (osklivost0>0.0)
   {
   	xprintf(Msg,"\nzeta0 se nepodarilo nastavit (osklivost je 1.0+(%f))", osklivost0-1.0);
   }
   free(zeta0);
   zeta0 = NULL;
   free(zeta);
   zeta = NULL;

   return vratit;
}

void che_maticovy_vypocet (char *soubor)
{
   int i = 0;
   int j = 0;
   int pruchodu = 0;
   double *zeta = NULL; //P_che[].zeta0;
   double *KK = NULL;//=K(rce);
   double *skala = NULL;
   double *J = NULL;
   double *hodnota = NULL;
   int *hprvky = NULL;
   double *delta = NULL;
   double *prstr = NULL;
   FILE *fw = NULL;
   double norma = 0.0;
   double stara_norma = 0.0;
   int stagnace = 0;

   if (G_prm.pocet_reakci_pro_matici == 0)
   {
      for (i=0; i<G_prm.pocet_latekvefazi; i++)
      {
         P_lat[i].m = P_lat[i].m0;
      }
      return;
   }
   if (che_urci_zeta0() > 0)
   {
      for (i=0; i<G_prm.pocet_latekvefazi; i++)
      {
         P_lat[i].m = P_lat[i].m0;
      }
      fw = fopen(soubor, "a");
      fprintf (fw,"\nchemie: URCITE NEKOREKTNI POCATECNI PODMINKA!\t");
      fclose(fw);
      return;
   }
   zeta = (double *)malloc( G_prm.pocet_reakci_pro_matici*sizeof( double ));
   if ( zeta == NULL )
   {
   	printf ("Malo pameti!\n");
   	exit(0);
   }for(j = 0; j < G_prm.pocet_reakci_pro_matici;j++){zeta[j] = 0.0;}

   prstr = (double *)malloc( G_prm.pocet_reakci_pro_matici*sizeof( double ));
   if ( prstr == NULL )
   {
   	printf ("Malo pameti!\n");
   	exit(0);
   }for(j = 0; j < G_prm.pocet_reakci_pro_matici;j++){prstr[j] = 0.0;}

   delta = (double *)malloc( G_prm.pocet_reakci_pro_matici*sizeof( double ));
   if ( delta == NULL )
   {
   	printf ("Malo pameti!\n");
   	exit(0);
   }for(j = 0; j < G_prm.pocet_reakci_pro_matici;j++){delta[j] = 0.0;}

   KK = (double *)malloc( G_prm.pocet_reakci_pro_matici*sizeof( double ));
   if ( KK == NULL )
   {
   	printf ("Malo pameti!\n");
   	exit(0);
   }for(j = 0; j < G_prm.pocet_reakci_pro_matici;j++){KK[j] = 0.0;}

   skala = (double *)malloc( G_prm.pocet_reakci_pro_matici*sizeof( double ));
   if ( skala == NULL )
   {
   	printf ("Malo pameti!\n");
   	exit(0);
   }for(j = 0; j < G_prm.pocet_reakci_pro_matici;j++){skala[j] = 0.0;}

   J = (double *)malloc( G_prm.pocet_reakci_pro_matici*G_prm.pocet_reakci_pro_matici*sizeof( double ));
   if ( J == NULL )
   {
   	printf ("Malo pameti!\n");
   	exit(0);
   }for(j = 0; j < G_prm.pocet_reakci_pro_matici;j++){J[j] = 0.0;}

   hodnota = (double *)malloc( G_prm.pocet_reakci_pro_matici*sizeof( double ));
   if ( hodnota == NULL )
   {
   	printf ("Malo pameti!\n");
   	exit(0);
   }for(j = 0; j < G_prm.pocet_reakci_pro_matici;j++){hodnota[j] = 0.0;}
   hprvky = (int *)malloc( G_prm.pocet_reakci_pro_matici*sizeof( int ));
   if ( hprvky == NULL )
   {
   	printf ("Malo pameti!\n");
   	exit(0);
   }for(j = 0; j < G_prm.pocet_reakci_pro_matici;j++){hprvky[j] = 0.0;}
   for ( i=0; i<G_prm.pocet_reakci_pro_matici; i++)
   {
      KK[i]=che_K_(i);
      zeta[i]=P_che[i].zeta0;
   }
   che_srovnej(KK,skala);
   if (G_prm.vypisy>1)
   {
	   fw = fopen(soubor, "a");
	   fprintf (fw,"\nK -> ");
	   che_napis_soubor_ld(fw, KK, G_prm.pocet_reakci_pro_matici);
	   fclose(fw);
   }

   che_hodnoty(zeta,hodnota,skala);
   che_Jakobi(zeta,J,skala);
   che_odecti(hodnota, KK, prstr);
   if (G_prm.vypisy>2)
   {
	   che_napis_stav_soubor(soubor,zeta, hodnota, J, prstr);
   }
   pruchodu = 0;
   stagnace = 0;
   stara_norma = che_norma(prstr,KK);
   while ((che_norma(prstr,KK)>G_prm.epsilon)&&(pruchodu < MAX_POC_VNEJ_CYK)&&(stagnace<MAX_STAGNACE)) 	// 2x cislo vycucane z prstu
   {
   	pruchodu++;
      norma = che_norma(prstr,KK);
//printf ("stara %Le nova %Le", stara_norma, norma);
      if (norma>stara_norma)
      {
         G_prm.omega /= 2.0;
         if (G_prm.vypisy>4)
         {
            fw = fopen(soubor, "a");
            fprintf(fw,"\n    ( - omega = %f )", G_prm.omega);
            fclose (fw);
         }
         stara_norma = norma;
         stagnace = 0;
      }
      else if (norma*1.5 < stara_norma) //cislo vycucane z prstu
      {
         G_prm.omega *= 2.0;
         if ( G_prm.omega > 1.5 )  //cislo vycucane z prstu
         {
            G_prm.omega = 1.5;      //cislo vycucane z prstu
         }
         if ( G_prm.omega < -1.5 )  //cislo vycucane z prstu
         {
            G_prm.omega = -1.5;     //cislo vycucane z prstu
         }
         if (G_prm.vypisy>4)
         {
            fw = fopen(soubor, "a");
            fprintf(fw,"\n    ( + omega = %f )", G_prm.omega);
            fclose (fw);
         }
         stara_norma = norma;
         stagnace = 0;
      }
      else
      {
         stagnace++;
      }
		if (G_prm.vypisy>4)
      {
		   fw = fopen(soubor, "a");
      	fprintf(fw,"\n%d. VNEJSI CYKLUS", pruchodu);
         fclose (fw);
   		if (G_prm.vypisy>5)
         {
            fw = fopen(soubor, "a");
            fprintf(fw,"\nmolality:");
            for (i=0;i<G_prm.pocet_latekvefazi;i++)
            {
               fprintf(fw,"\t%f",che_m(i,zeta));
            }
            fclose (fw);
         }
      }
	   if (che_Gauss ( J, prstr, hprvky, G_prm.pocet_reakci_pro_matici) > 0)
	   {
		   fw = fopen(soubor, "a");
		   fprintf (fw,"\nchemie: Gauss spadnul!");
  			fclose(fw);
  		   exit(222);
	   }
	   che_zgaussproprg (prstr,hprvky,delta);
		if (G_prm.vypisy>2)
      {
		   fw = fopen(soubor, "a");
		   fprintf (fw, "\n>  delta=");
		   che_napis_soubor_ld(fw, delta, G_prm.pocet_reakci_pro_matici);
		   fclose(fw);
      }
	   che_nasob_ld(G_prm.omega,delta,delta, G_prm.pocet_reakci_pro_matici);
	   if (che_odecti_s_korekci_ld (zeta,delta,zeta,G_prm.pocet_reakci_pro_matici) > 0)
      {
		   fw = fopen(soubor, "a");
			fprintf (fw,"\nchemie: PATRNE NEKOREKTNI POCATECNI PODMINKA!\t");
  			fclose(fw);
         che_nuluj_ld(zeta, G_prm.pocet_reakci_pro_matici);
      	break;
      }
	   che_hodnoty(zeta,hodnota,skala);
	   che_Jakobi(zeta,J,skala);
	   che_odecti(hodnota, KK, prstr);
	   if (G_prm.vypisy>2)
      {
		   che_napis_stav_soubor(soubor, zeta, hodnota, J, prstr);
      }
   }
   if (stagnace >= MAX_STAGNACE)					 						   // cislo vycucane z prstu
   {
	   fw = fopen(soubor, "a");
      fprintf (fw, "\nchemie: NEMOHU DODRZET POZADOVANOU PRESNOST!\n");
		fclose(fw);
   }
   if (pruchodu >= MAX_POC_VNEJ_CYK)					 						   // cislo vycucane z prstu
   {
	   fw = fopen(soubor, "a");
      fprintf (fw, "\nchemie: PATRNE PRILIS RYCHLE KINETICKE REAKCE!\t");
		fclose(fw);
	   exit(223);
   }
   for (i=0; i<G_prm.pocet_latekvefazi; i++)
   {
      P_lat[i].m = che_m(i,zeta);
      P_lat[i].m0 = che_m(i,zeta);
   }
   free(KK);
   KK = NULL;
   free(zeta);
   zeta = NULL;
   free(delta);
   delta = NULL;
   free(prstr);
   prstr = NULL;
   free(hprvky);
   hprvky = NULL;
   free(hodnota);
   hodnota = NULL;
   free(J);
   J = NULL;
   free(skala);
   skala = NULL;
}

void che_spocitej_rychlosti(FILE *fw, double *rychlosti, double *poloha, double dt)
{
	int i = 0;
	int j = 0;
   int chyba = 0;

   for (i=G_prm.pocet_reakci_pro_matici+G_prm.pocet_rozpadu;i<G_prm.celkovy_pocet_reakci;i++)
   {
//     if (G_prm.vypisy > 4)
//       {
// 	      fprintf(fw,"\n    rychlost %d. kineticke reakce a casovy krok", i);
//       }
   	rychlosti[i-G_prm.pocet_reakci_pro_matici-G_prm.pocet_rozpadu] = P_che[i].K*dt;
      for (j=0;j < G_prm.pocet_latekvefazi;j++)
      {
	   	rychlosti[i-G_prm.pocet_reakci_pro_matici-G_prm.pocet_rozpadu] *= che_poww_ld(poloha[j],P_che[i].exponent[j],&chyba);
	 if (G_prm.vypisy > 4)
	 {
	      fprintf(fw,"\n  %d.  rychlost %d. kineticke reakce %10.24f, poloha = %f, exponent = %f\n", j, i, rychlosti[i-G_prm.pocet_reakci_pro_matici-G_prm.pocet_rozpadu], poloha[j], P_che[i].exponent[j]);
	 }
         if (chyba > 0)
         {
         	// nula na zaporne cislo nebo zaporne cislo na necele cislo!
         }
      }
		if (G_prm.vypisy > 4)
      {
//       	fprintf(fw," = %f. %f", rychlosti[i-G_prm.pocet_reakci_pro_matici-G_prm.pocet_rozpadu],dt);
		}
   }
}

void che_spocitej_posunuti (double *posunuti, double *rychlosti)
{
	int i = 0;
	int j = 0;

   che_nuluj_ld(posunuti,G_prm.pocet_latekvefazi);
   for (i=G_prm.pocet_reakci_pro_matici+G_prm.pocet_rozpadu;i<G_prm.celkovy_pocet_reakci;i++)
   {
      for (j=0;j<G_prm.pocet_latekvefazi;j++)
      {
         posunuti[j] += P_che[i].stech_koef_p[j]*rychlosti[i-G_prm.pocet_reakci_pro_matici-G_prm.pocet_rozpadu];
      }
	}
}

void che_prepocitej_polohu (double *poloha2, double *poloha, double *posunuti)
{
   int j = 0;

   for (j=0;j<G_prm.pocet_latekvefazi;j++)
   {
      poloha2[j] = poloha[j]+posunuti[j];
   }
}

void che_zkrat_latku_o(FILE *fw, int kterou, double o_kolik, double *rychlosti)
{
// mozna bych mel vsechny reakce spotrebovavajici latku zkratit rovnym dilem.
   int i = 0;
   int reakce = 0;
   double maximum = 0.0;

   for (;;)
   {
      reakce = -1;
      maximum = 0.0;
      for (i=G_prm.pocet_reakci_pro_matici+G_prm.pocet_rozpadu;i<G_prm.celkovy_pocet_reakci;i++)
      {
         if (-P_che[i].stech_koef_p[kterou]*rychlosti[i-G_prm.pocet_reakci_pro_matici-G_prm.pocet_rozpadu]>maximum)
         {
            maximum = -P_che[i].stech_koef_p[kterou]*rychlosti[i-G_prm.pocet_reakci_pro_matici-G_prm.pocet_rozpadu];
            reakce = i;
         }
      }
      if (reakce==-1)
      {
         fprintf (fw,"\nchemie: Tohle se vubec nemelo stat, nerozumim tomu - nelze uz zkratit spotrebu latky!");
         fclose(fw);
         exit(224);
      }
      if (maximum>=o_kolik)
      {
         rychlosti[reakce-G_prm.pocet_reakci_pro_matici-G_prm.pocet_rozpadu] += o_kolik/P_che[reakce].stech_koef_p[kterou];
         return;
      }
      else
      {
         o_kolik -= maximum;
         rychlosti[reakce-G_prm.pocet_reakci_pro_matici-G_prm.pocet_rozpadu] = 0.0;
      }
   }
}

void che_pomala_kinetika (char *soubor, int poc_kroku)
{
   FILE *fw = NULL;
   int i,j, krok;
   double *poloha = NULL;
   double *poloha2 = NULL;
   double *k1 = NULL;
   double *k2 = NULL;
   double *k3 = NULL;
   double *k4 = NULL;
   double *rychlosti = NULL;
   double *posunuti = NULL;
   double dt = 0.0;

   if (G_prm.pocet_reakci_pro_matici+G_prm.pocet_rozpadu+G_prm.pocet_pom_kin != G_prm.celkovy_pocet_reakci)
   {
      printf ("POMALA KINETIKA POZADUJE, ABY NEEXISTOVALY JINE REAKCE NEZ PRO MATICI, POMALE KINETICKE A ROZPADY! %d %d %d %d !", G_prm.pocet_reakci_pro_matici,G_prm.pocet_rozpadu,G_prm.pocet_pom_kin,G_prm.celkovy_pocet_reakci);
      exit (115);
   }
   if (G_prm.pocet_pom_kin==0)
   {
      return;
   }
   fw = fopen(soubor, "a");
   if (G_prm.vypisy>4) xprintf(Msg,"\nche_pomala_kinetika: ");
/*
//  TOHLE JE "ZPARCHANTELEJ EULER" - EULER PRO KAZDOU ROVNICI ZVLAST
	long int i,j;
   double rychlost;

  	for (i=G_prm.pocet_reakci_pro_matici;i<G_prm.celkovy_pocet_reakci;i++)
   {
		if (G_prm.vypisy>4)
      {
         fprintf (fw,"\n    rychlost %ld. kineticke reakce", i);
      }
   	rychlost = P_che[i].K*G_prm.deltaT;
      for (j=0;j<G_prm.pocet_latekvefazi;j++)
      {
         rychlost *= che_poww_ld(P_lat[j].m0,P_che[i].exponent[j]);
      }
		if (G_prm.vypisy>4)
      {
      	fprintf (fw," = %Le.", rychlost);
		}
      for (j=0;j<G_prm.pocet_latekvefazi;j++)
      {
      	P_lat[j].m0 += P_che[i].stech_koef_p[j]*rychlost;
         if (P_lat[j].m0<0.0)
         {
            fprintf (fw,"\nchemie: Pomale kineticke reakce nejsou dost pomale!");
            fclose(fw);
            exit(224);
         }
      }
   }
*/
/*
//  TOHLE JE EULER
	long int i,j;
	long double *rychlosti;
   long double *posunuti;

	rychlosti = (long double *)malloc( (G_prm.celkovy_pocet_reakci-G_prm.pocet_reakci_pro_matici)*sizeof( long double ));
   if ( rychlosti == NULL )
   {
   	printf ("Malo pameti!\n");
   	exit(0);
   }
   for (i=G_prm.pocet_reakci_pro_matici;i<G_prm.celkovy_pocet_reakci;i++)
   {
		if (G_prm.vypisy>4)
      {
	      fprintf (fw,"\n    rychlost %ld. kineticke reakce", i);
      }
   	rychlosti[i-G_prm.pocet_reakci_pro_matici] = P_che[i].K*G_prm.deltaT;
      for (j=0;j<G_prm.pocet_latekvefazi;j++)
      {
	   	rychlosti[i-G_prm.pocet_reakci_pro_matici] *= che_poww_ld(P_lat[j].m0,P_che[i].exponent[j]);
      }
		if (G_prm.vypisy>4)
      {
      	fprintf (fw," = %Le.", rychlosti[i-G_prm.pocet_reakci_pro_matici]);
		}
   }
	posunuti = (long double *)malloc( G_prm.pocet_latekvefazi*sizeof( long double ));
   if ( posunuti == NULL )
   {
   	printf ("Malo pameti!\n");
   	exit(0);
   }
   che_nuluj_ld(posunuti,G_prm.pocet_latekvefazi);
   for (i=G_prm.pocet_reakci_pro_matici;i<G_prm.celkovy_pocet_reakci;i++)
   {
      for (j=0;j<G_prm.pocet_latekvefazi;j++)
      {
         posunuti[j] += P_che[i].stech_koef_p[j]*rychlosti[i-G_prm.pocet_reakci_pro_matici];
      }
	}
   for (j=0;j<G_prm.pocet_latekvefazi;j++)
   {
      P_lat[j].m0 += posunuti[j];
      if (P_lat[j].m0<0.0)
      {
         fprintf (fw,"\nchemie: Pomale kineticke reakce nejsou dost pomale!");
         fclose(fw);
         exit(224);
      }
   }
   free(rychlosti);
   free(posunuti);
*/
//  TOHLE JE JEDEN KROK RUNGE-KUTTA - mel bych priprogramovat moznost rozdelit vypocet na vic kroku

	poloha = (double *)malloc( (G_prm.pocet_latekvefazi)*sizeof( double ));
   if ( poloha == NULL )
   {
   	printf ("Malo pameti!\n");
   	exit(0);
   }for(i = 0;i < G_prm.pocet_latekvefazi; i++){ poloha[i] = 0.0; }

	poloha2 = (double *)malloc( (G_prm.pocet_latekvefazi)*sizeof( double ));
   if ( poloha2 == NULL )
   {
   	printf ("Malo pameti!\n");
   	exit(0);
   }for(i = 0;i < G_prm.pocet_latekvefazi; i++){ poloha2[i] = 0.0; }

	rychlosti = (double *)malloc( G_prm.pocet_pom_kin*sizeof( double ));
   if ( rychlosti == NULL )
   {
   	printf ("Malo pameti!\n");
   	exit(0);
   }for(i = 0;i < G_prm.pocet_pom_kin; i++){ rychlosti[i] = 0.0; }

	k1 = (double *)malloc( G_prm.pocet_pom_kin*sizeof( double ));
   if ( k1 == NULL )
   {
   	printf ("Malo pameti!\n");
   	exit(0);
   }for(i = 0;i < G_prm.pocet_pom_kin; i++){ k1[i] = 0.0; }

	k2 = (double *)malloc( G_prm.pocet_pom_kin*sizeof( double ));
   if ( k2 == NULL )
   {
   	printf ("Malo pameti!\n");
   	exit(0);
   }for(i = 0;i < G_prm.pocet_pom_kin; i++){ k2[i] = 0.0; }

	k3 = (double *)malloc( G_prm.pocet_pom_kin*sizeof( double ));
   if ( k3 == NULL )
   {
   	printf ("Malo pameti!\n");
   	exit(0);
   }for(i = 0;i < G_prm.pocet_pom_kin; i++){ k3[i] = 0.0; }

	k4 = (double *)malloc( G_prm.pocet_pom_kin*sizeof( double ));
   if ( k4 == NULL )
   {
   	printf ("Malo pameti!\n");
   	exit(0);
   }for(i = 0;i < G_prm.pocet_pom_kin; i++){ k4[i] = 0.0; }

	posunuti = (double *)malloc( G_prm.pocet_latekvefazi*sizeof( double ));
   if ( posunuti == NULL )
   {
   	printf ("Malo pameti!\n");
   	exit(0);
   }for(i = 0;i < G_prm.pocet_latekvefazi; i++){ posunuti[i] = 0.0; }
   dt= G_prm.deltaT / poc_kroku;
   for (krok = 0; krok<poc_kroku; krok++)
   {
      for (j=0;j<G_prm.pocet_latekvefazi;j++)
      {
         poloha[j] = P_lat[j].m0;
         if (P_lat[j].m0<0.0)
         {
            fprintf (fw,"\nchemie: Vstup do pomalych kinetickych reakci obsahoval zapornou molalitu %d. latky-neprobehne vypocet!!", j+1);
            free(rychlosti);
	    rychlosti = NULL;
            
	    free(k1);
	    k1 = NULL;

            free(k2);
	    k2 = NULL;

            free(k3);
	    k3 = NULL;

            free(k4);
	    k4 = NULL;

            free(posunuti);
	    posunuti = NULL;

            free(poloha2);
	    poloha2 = NULL;

            free(poloha);
	    poloha = NULL;
            return;
         }
      }
      che_spocitej_rychlosti(fw, k1, poloha, dt);
      che_nasob_ld (0.5, k1, rychlosti, G_prm.pocet_pom_kin);
      che_spocitej_posunuti(posunuti, rychlosti);
      che_prepocitej_polohu(poloha2, poloha, posunuti);

      che_spocitej_rychlosti(fw, k2, poloha2, dt);
      che_nasob_ld (0.5, k2, rychlosti, G_prm.pocet_pom_kin);
      che_spocitej_posunuti(posunuti, rychlosti);
      che_prepocitej_polohu(poloha2, poloha, posunuti);

      che_spocitej_rychlosti(fw, k3, poloha2, dt);
      che_spocitej_posunuti(posunuti, k3);
      che_prepocitej_polohu(poloha2, poloha, posunuti);

      che_spocitej_rychlosti(fw, k4, poloha2, dt);

      che_kombinuj4_ld(1.0/6.0, k1, 1.0/3.0, k2, 1.0/3.0, k3, 1.0/6.0, k4, rychlosti, G_prm.pocet_pom_kin);
      che_spocitej_posunuti(posunuti, rychlosti);
      che_prepocitej_polohu(poloha2, poloha, posunuti);

      for (j=0;j<G_prm.pocet_latekvefazi;j++)
      {
         if (poloha2[j]<0.0)
         {
            fprintf (fw,"\nchemie: pri pomalych kinetickych reakcich dosla %d. latka (%f)\t", j+1, poloha2[j]);
   //         fclose(fw);
   //         exit(224);
            che_zkrat_latku_o(fw,j,-poloha2[j],rychlosti);
            che_spocitej_posunuti(posunuti, rychlosti);
            che_prepocitej_polohu(poloha2, poloha, posunuti);
         }
      }
      for (j=0;j<G_prm.pocet_latekvefazi;j++)
      {
         P_lat[j].m0 = poloha2[j];
         P_lat[j].m = poloha2[j];
         if (P_lat[j].m0<0.0)
         {
            if (P_lat[j].m0>-1.e-20)						// cislo vycucane z prstu
            {
               P_lat[j].m0 = 0.0;
            }
            else
            {
               fprintf (fw,"\nchemie: Tohle se vubec nemelo stat, nerozumim tomu - pomale kineticke reakce nejsou dost pomale!\n%d.latka (%f)", j+1,P_lat[j].m0);
               fclose(fw);
               exit(224);
            }
         }
      }
   }
   free(rychlosti);
   rychlosti = NULL;
            
   free(k1);
   k1 = NULL;

   free(k2);
   k2 = NULL;

   free(k3);
   k3 = NULL;

   free(k4);
   k4 = NULL;

   free(posunuti);
   posunuti = NULL;

   free(poloha2);
   poloha2 = NULL;

  free(poloha);
  poloha = NULL;
//
   fclose (fw);
  	if (G_prm.vypisy>4) xprintf(Msg,"o.k. (che_pomala_kinetika)");
}

void che_rovnovazne_sorpce (char *soubor)
{
   FILE *fw = NULL;
   int j = 0;
   double celk_lat_mnoz = 0.0;
   int chyba = 0;
   double mmin = 0.0;
   double mmax = 0.0;
   double mpul = 0.0;
   double spul = 0.0;
   double xmin = 0.0; 
   double xmax = 0.0;
   double xpul = 0.0;

   fw = fopen(soubor, "a");
   if (G_prm.vypisy>4) xprintf(Msg,"\nche_sorpce: ");
   for (j=0; j<G_prm.pocet_latekvefazi; j++)
   {
      P_lat[j].m_sorb = P_lat[j].m0_sorb;
//printf("\nTyp sorpce je %ld",P_lat[j].typ_sorpce);
      if (P_lat[j].typ_sorpce == 0) continue;
      celk_lat_mnoz=P_lat[j].m0_sorb*G_prm.splocha+P_lat[j].m*G_prm.objem;
      if (celk_lat_mnoz > 1e-16)
      {
         switch (P_lat[j].typ_sorpce)
         {
            case 1: //linearni
//printf("\nParametr linearni izotermy je %Lf",P_lat[j].param_sorpce[0]);
               P_lat[j].m=celk_lat_mnoz/(G_prm.objem+G_prm.splocha*P_lat[j].param_sorpce[0]);
               break;
            case 2: //Freundlichova
//printf("\nParametry Freundlichovy izotermy jsou %Lf, %Lf",P_lat[j].param_sorpce[0],P_lat[j].param_sorpce[1]);
               {mmin = 0.0;
               xmin = 0.0;
               mmax = celk_lat_mnoz/G_prm.objem;
			   xmax = celk_lat_mnoz + G_prm.splocha*P_lat[j].param_sorpce[0]*che_poww_ld(mmax,P_lat[j].param_sorpce[1],&chyba);
               if (chyba > 0)
               {
                  xprintf(Msg,"\nZcela nepochopitelne nemuzu mocnit pri vypoctu Freundlichovy izotermy!");
                  exit(100);
               }
               do{
//printf("\n(x= %Le, xmin = %Le, xmax = %Le)", celk_lat_mnoz, xmin, xmax);
                  mpul = mmin+(celk_lat_mnoz-xmin)/(xmax-xmin)*(mmax-mmin);
				  spul = P_lat[j].param_sorpce[0]*che_poww_ld(mpul,P_lat[j].param_sorpce[1],&chyba);
                  if (chyba > 0)
                  {
                     xprintf(Msg,"\nZcela nepochopitelne nemuzu mocnit pri vypoctu Freundlichovy izotermy!");
                     exit(100);
                  }
                  xpul = mpul * G_prm.objem + spul * G_prm.splocha;
                  if (xpul>celk_lat_mnoz)
                  {
                     mmax = mpul;
                     xmax = xpul;
                  }
                  else
                  {
                     mmin = mpul;
                     xmin = xpul;
                  }
//printf(" mmin = %Le, mmax = %Le, mpul = %Le, spul = %Le", mmin, mmax, mpul, spul);
               } while (fabs(xpul-celk_lat_mnoz)>celk_lat_mnoz*1.0e-8);      //cislo vycucane z prstu
               P_lat[j].m = mpul;}
               break;
            case 3: //Langmuirova
//printf("\nParametry Langmuirovy izotermy jsou %Lf, %Lf",P_lat[j].param_sorpce[0],P_lat[j].param_sorpce[1]);
               P_lat[j].m=(P_lat[j].param_sorpce[0]*celk_lat_mnoz-G_prm.objem-G_prm.splocha*P_lat[j].param_sorpce[0]*P_lat[j].param_sorpce[1]
						  +sqrt(che_poww(G_prm.splocha*P_lat[j].param_sorpce[0]*P_lat[j].param_sorpce[1]+G_prm.objem-P_lat[j].param_sorpce[0]*celk_lat_mnoz,2,&chyba)+4.0*celk_lat_mnoz*G_prm.objem*P_lat[j].param_sorpce[0]))
                          /2.0/G_prm.objem/P_lat[j].param_sorpce[0];
//printf(" (celk_lat_mnoz= %Le, m = %Le)", celk_lat_mnoz, P_lat[j].m);
               if (chyba > 0)
               {
                  xprintf(Msg,"\nZcela nepochopitelne nemuzu mocnit pri vypoctu Langmuirovy izotermy!");
                  exit(100);
               }
               break;
            default: //neznama
               xprintf(Msg,"\nNeznamy typ sorpce (cislo %d)!", P_lat[j].typ_sorpce);
               exit(100);
         }
   	}
      else
      {
	      P_lat[j].m = 0.0;
      }
      P_lat[j].m_sorb = (celk_lat_mnoz-P_lat[j].m*G_prm.objem)/G_prm.splocha;
      P_lat[j].m0_sorb = P_lat[j].m_sorb;
      P_lat[j].m0 = P_lat[j].m;
   }
   fclose (fw);
  	if (G_prm.vypisy>4) xprintf(Msg,"o.k. (che_sorpce)");
}

void che_rozpad (char *soubor)
{
   FILE *fw = NULL;
   int i = 0;
   int j = 0;
   double dt = 0.0;
   double pomer = 0.0;
   double m0 = 0.0;
   double rychlost = 0.0;

   if (G_prm.pocet_rozpadu==0)
   {
      return;
   }
   fw = fopen(soubor, "a");
   if (G_prm.vypisy>4) xprintf(Msg,"\nche_rozpad: ");
   dt= G_prm.deltaT;
   for (j=0;j<G_prm.pocet_latekvefazi;j++)
   {
      P_lat[j].m = P_lat[j].m0;
   }
   for (i=G_prm.pocet_reakci_pro_matici; i<G_prm.pocet_reakci_pro_matici+G_prm.pocet_rozpadu; i++)
   {
      pomer = 1.0-exp(-log(2.0)/P_che[i].K*dt);
//printf("\npomer = 1.0-expl(%Lf/%Le*%Le=%Le)",-logl(2.0),P_che[i].K,dt,pomer);
      m0=-10.0;
      for (j=0;j<G_prm.pocet_latekvefazi;j++)
      {
         if (P_che[i].stech_koef_p[j]==-1)
         {
            m0 = P_lat[j].m0;
            break;
         }
      }
      if (m0 == -10.0)
      {
         printf ("\nV nektere rozpadove rakci se nerozpada zadna latka se stech. koef. -1!");
         exit (133);
      }
      rychlost = pomer*m0;
      for (j=0;j<G_prm.pocet_latekvefazi;j++)
      {
         P_lat[j].m += P_che[i].stech_koef_p[j]*rychlost;
      }
   }
   for (j=0;j<G_prm.pocet_latekvefazi;j++)
   {
      P_lat[j].m0 = P_lat[j].m;
   }
   fclose (fw);
  	if (G_prm.vypisy>4) xprintf(Msg,"o.k. (che_rozpad)");
}

void che_kineticke_sorpce (char *soubor)
{
////////////////////// SEM VRAZIT KINETICKE SORPCE ///////////////////////
//                                                                      //
//                                                                      //
//                                                                      //
//                                                                      //
//                                                                      //
//                                                                      //
//                                                                      //
////////////////////// SEM VRAZIT KINETICKE SORPCE ///////////////////////
}

void che_vypocet_rovnovah (char *soubor)
{
   G_prm.pocet_reakci_pro_matici = G_prm.pocet_rovnovah;
   che_maticovy_vypocet(soubor);
   G_prm.pocet_reakci_pro_matici = G_prm.pocet_rovnovah+G_prm.pocet_kinetik;
}

void che_rovnovahy_se_sorpci(char *soubor)
{
// sem muzu vrazit cyklus
   che_vypocet_rovnovah(soubor);
   che_rovnovazne_sorpce(soubor);
   che_presun_poc_p_();
}

void che_matice_se_sorpci(char *soubor)
{
// sem nesmim vrazit cyklus - jsou tam i kineticke reakce
   che_maticovy_vypocet(soubor);
   che_rovnovazne_sorpce(soubor);
   che_presun_poc_p_();
}

void che_pocitej_soubor(char *soubor, int *poc_krok)
{
	if (poc_krok > 0)
   {
//// tohle zaremovani doufam pomuze na rychle kineticke reakce v matici
//   	if (G_prm.celkovy_pocet_reakci>G_prm.pocet_reakci_pro_matici)
      {
      	che_rovnovahy_se_sorpci(soubor);
      }
   	poc_krok = 0;
//      che_vypis_soubor(soubor);
   }
	che_pomala_kinetika(soubor,G_prm.deleni_RK);
	che_rozpad(soubor);
   //che_kineticke_sorpce(soubor);
// sem muzu pridat che_rovnovahy_se_sorpci(soubor);
	che_matice_se_sorpci(soubor);
}

void che_nadpis__soubor(char *soubor)
{
	int i = 0;
   FILE *fw = NULL;

   fw = fopen(soubor, "a");
	fprintf (fw,"\nkrok\tcas");
   for (i=0; i<G_prm.pocet_latekvefazi; i++)
   {
 		fprintf (fw,"\t%d. latka(rozp.)", i);
   }
   for (i=0; i<G_prm.pocet_latekvefazi; i++)
   {
 		fprintf (fw,"\t%d. latka(sorb.)", i);
   }
   fprintf(fw,"\tobjem");
	fprintf (fw,"\n0\t0.0");
   che_outpocp__soubor(fw);
   fclose(fw);
}

//*************************************************************************
//                 FUNKCE NA NACTENI CHEMIE PRO TRANSPORT
//*************************************************************************
void che_nactenichemie(void)
{
}

void che_vypis_prm_lat_che ( void )
{
   int i = 0;
   int j = 0;
   xprintf(Msg,"\nPRM: ");
   xprintf(Msg,"%s %d %d %d %d ", G_prm.jmeno_ich, G_prm.pocet_latek,G_prm.pocet_latekvefazi,G_prm.celkovy_pocet_reakci,G_prm.pocet_reakci_pro_matici);
   xprintf(Msg,"%f %f %f %f %f %f ",G_prm.T,G_prm.Afi,G_prm.b,G_prm.epsilon,G_prm.omega,G_prm.deltaT);
   xprintf(Msg,"%d %d",G_prm.cas_kroku,G_prm.vypisy);

   xprintf(Msg,"\nLAT: ");
   for (i=0; i<G_prm.pocet_latek; i++)
   {
       xprintf(Msg,"\n  (%d): ", i);
//        xprintf(Msg,"%s %f %f %f %f %d %f",P_lat[i].nazev,P_lat[i].m0,P_lat[i].m,P_lat[i].M,P_lat[i].dGf,P_lat[i].Q,P_lat[i].aktivita);
       xprintf(Msg,"%d. %f %f %f %f %d %f",i,P_lat[i].m0,P_lat[i].m,P_lat[i].M,P_lat[i].dGf,P_lat[i].Q,P_lat[i].aktivita);
   }

   xprintf(Msg,"\nCHE: ");
   for (i=0; i<G_prm.celkovy_pocet_reakci; i++)
   {
       xprintf(Msg,"\n  (%d): ", i);
//        xprintf(Msg,"%s %f %d %f",P_che[i].nazev,P_che[i].K,P_che[i].typ_reakce,P_che[i].zeta0);
       xprintf(Msg,"%d %f %d %f",i,P_che[i].K,P_che[i].typ_reakce,P_che[i].zeta0);
       xprintf(Msg,"\n     stech. koef.: ");
      for (j = 0; j<G_prm.pocet_latek; j++)
      {
          xprintf(Msg,"%d ", P_che[i].stech_koef_p[j]);
      }
       xprintf(Msg,"\n     exponenty:    ");
      for (j = 0; j<G_prm.pocet_latekvefazi; j++)
      {
          xprintf(Msg,"%f ", P_che[i].exponent[j]);
      }
   }
}

/********************************************************************/
/*						Uvolnuje pamet     									     */
/********************************************************************/
void che_uvolneni_P( void )
{
   xprintf(Msg, "Uvolneni che_ P_lat, P_che : " );
   if ( P_lat != NULL )
   {
      free( P_lat ); P_lat = NULL;
   }
   xprintf(Msg, "O.k., " );
   if ( P_che != NULL )
   {
      free( P_che ); P_che = NULL;
   }
   xprintf(Msg, "O.k.\n" );
}
