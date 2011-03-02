//---------------------------------------------------------------------------
#ifndef interfaceH
#define interfaceH
//---------------------------------------------------------------------------
#ifndef mobile
#define mobile 0
#endif
//---------------------------------------------------------------------------
#ifndef immobile
#define immobile 1
#endif
//---------------------------------------------------------------------------
//#ifndef mobile_sorb
//#define mobile_sorb 2
//#endif
//---------------------------------------------------------------------------
//#ifndef immobile_sorb
//#define immobile_sorb 3
//#endif
//---------------------------------------------------------------------------

#include <stdio.h>
#include <stdio.h>
#include <dirent.h>
#include "che_head.h"
void print_usage(void);
void kontrola(void);
void priprav(void);

//--------------pro semchem--------------------------------------------------
void che_nadpis__soubor(char *soubor);
void che_outpocp_soubor(FILE *fw);
void che_pocitej_soubor(char *soubor, int *poc_krok);
void che_vypis_soubor(char *soubor);
void che_presun_poc_p_(void);
void che_vypis__soubor(char *soubor);

//--------------------------------------------------------------------------
//  GLOBALNI PROMENNE
//---------------------------------------------------------------------------
extern struct TS_prm	G_prm;
extern struct TS_lat 	*P_lat;
extern struct TS_che	*P_che;

#endif
