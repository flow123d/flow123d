#include <stdlib.h>                                                   // ! Dressler 24.08.2001
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <dirent.h>
#include <petscvec.h>
#include "elements.h"

#ifdef MAIN
#define EX
#else
#define EX extern
#endif

#ifndef MAXPATH
#define MAXPATH 2048
#endif

#define MAX_POC_LATEK 200
#define MAX_POC_PARAM_SORPCE 5

#define ACCURACY 1e-12

struct TS_prm //typedef struct
{
   char jmeno_ich[ MAXPATH ];
	int	pocet_latek;
   int pocet_latekvefazi;
	int	celkovy_pocet_reakci;
	int	pocet_reakci_pro_matici;
	int	pocet_rovnovah;
	int	pocet_kinetik;
	int	pocet_pom_kin;
	int	pocet_rozpadu;
     double T; //long long double
     double TGf; //long long double
     double Afi; //long double
     double b;  //long double
      double epsilon; //long double
      double omega; //long double
      double deltaT;  //long double
      double objem; //long double
      double splocha; //long double
   int cas_kroku;
   int vypisy;
   int deleni_RK;
   int skaluj_matici;
   int abs_norma;
}; //TS_prm;

struct TS_lat //typedef struct
{
	  char nazev[ 80 ];      /* Jmeno latky */
	  double m0;        /* pocatecni molalita *///long double
//	  double m0_sorb;   /* pocatecni obsah sorbovany *///long double
	   double m;         /* konecna molalita *///long double
//	   double m_sorb;    /* konecny obsah sorbovany *///long double
	   double M;	        /* molarni hmotnost *///long double
	   double dGf;       /* prispevek ke Gibbsove energii *///long double
	   double dHf;       /* prispevek k entalpii *///long double
	int Q;            /* naboj */
	  double aktivita;  /* aktivita *///long double
	int typ_sorpce;   /* typ sorpce */
	   double param_sorpce[MAX_POC_PARAM_SORPCE];//long double
}; //TS_lat;

struct TS_che //typedef struct
{
	char nazev[ 80 ];		/* Jmeno reakce              */
	int stech_koef_p[MAX_POC_LATEK];	/* Stechiom.koef.slozek vpravo*/
	   double K;			/*rovnovazna, nebo kineticka konstanta, nebo bilance*///long double
	int typ_reakce;		/* 0==rovnovazna, 1==kineticka, 2==bilancni, 3==pomala kinetika */
	   double exponent[MAX_POC_LATEK];		/* exponenty pro kinetiku*///long double
      double zeta0;     /*pocatecni posunuti reakce*///long double
}; //TS_che;

//---------------------------------------------------------------------------
//  Funkce z che_semchem.cpp
//---------------------------------------------------------------------------
void che_nactenichemie( void );	     /* funkce nacteni chemie */
void ctiich (void); /*pomocna funkce nacteni chemie*/
//void che_vypocetchemie(struct Problem *problem, double **conc_mob_arr, double **conc_immob_arr, double **sorb_mob_arr, double **sorb_immob_arr);
void che_vypocetchemie(bool porTyp, double time_step, ElementIter ppelm, int poradi, double **conc_mob_arr, double **conc_immob_arr);
void che_uvolneni_P( void ); /* funkce na uvolneni pameti */

