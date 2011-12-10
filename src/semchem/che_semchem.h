#ifndef BIGHEAD
#define BIGHEAD
/********************************************************************/
/*   BIG-HEAD.H     GEN-TRAN Verse 1.0                              */
/********************************************************************/
#define DOS                                 /* Prepinac DOS X UNIX  */
#define VERB                                /* Upovidany vystup     */
#define BOOL int

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

#define fnsplit _splitpath                                            // ! Dressler 24.08.2001
#define fnmerge _makepath                                             // ! Dressler 24.08.2001
#ifdef DOS
 #include <dirent.h>
#endif
#include <stdlib.h>                                                   // ! Dressler 24.08.2001
#include <math.h>
#include <stdio.h>
#include <time.h>
#ifdef MAIN
 #define EX
#else
 #define EX extern
#endif

/********************************************************************/

// mh: kompiluju flow/tran
#define PROGRAM_TRAN

/********************************************************************/

#ifndef SKUPINY_MATERIALU
	#define SKUPINY_MATERIALU
#endif //SKUPINY_MATERIALU
/*------------------------------------------------------------------*/

#ifndef SKUPINA_SCENAR
	#define SKUPINA_SCENAR
#endif // SKUPINA_SCENAR
/*------------------------------------------------------------------*/

#ifndef SKUPINA_STENY
 #define SKUPINA_STENY
#endif // SKUPINA_STENY
/*------------------------------------------------------------------*/

#ifndef SKUPINA_RESICE
	#define SKUPINA_RESICE
#endif // SKUPINA_RESICE
/*------------------------------------------------------------------*/

#ifndef SKUPINA_TRANSPORT
 #define SKUPINA_TRANSPORT
#endif // SKUPINA_TRANSPORT
/*------------------------------------------------------------------*/

#ifndef SKUPINA_DUALPOROSITY
	#define SKUPINA_DUALPOROSITY
#endif // SKUPINA_DUALPOROSITY

/*------------------------------------------------------------------*/

#ifndef SKUPINA_REAKCE
	#define SKUPINA_REAKCE
#endif // SKUPINA_REAKCE
/********************************************************************/
#define MAXK1LM1         2      /* Max. pocet koeficientu v MELM    */
#define MAXELKOEF       10      /* Max. pocet koeficientu v ELM     */
#define MAXMATRKOEF      7      /* Max. pocet koef. v matr.         */
#define MAXDPORKOEF      3      /* Max. pocet koef. v DUAL POROSITY */
/********************************************************************/
#define NLENSPO         16      /* Pocet znaku popisu slozky R + H  */
#define NLENPVR          5      /* Maximalni delka popisu vrstev    */
/********************************************************************/
#define RUN_OK           0
#define RUN_ERROR        1
/********************************************************************/
// mh: kompiluju flow/tran
#define PRG_NAME         "tran"
/********************************************************************/
/********************************************************************/
/*                      Definice globalnich promennych              */
/********************************************************************/
/*  Vseobecne paramatry                                             */
/********************************************************************/
EX int      G_argc;
EX char     **G_argv;
EX char     *G_Program_name;
EX time_t   G_start_time;
EX char     *G_mezi_cas;
EX int      NULA;
EX int      G_exit_code;
/*------------------------------------------------------------------*/
#ifdef SKUPINA_RESICE
	EX int      G_cas_resice, G_cas_rozhrani;
#endif // SKUPINA_RESICE
/*------------------------------------------------------------------*/
#ifdef SKUPINA_CHEMIE
EX int      G_cas_chemie;
#endif // SKUPINA_CHEMIE
/********************************************************************/
/*   Seznam sousedu danych sten      ( -1 = vnejsi stena )          */
/********************************************************************/
#ifdef SKUPINA_STENY
	EX int      *P_Steny;
#endif // SKUPINA_STENY
/********************************************************************/
/*  Vektory pro ulozeni soustavy linearnich rovnic                  */
/********************************************************************/
#ifdef SKUPINA_RESICE
	EX int      *P_i, *P_j;
	EX double   *P_MM, *P_MR, *P_MX;
#endif // SKUPINA_RESICE
/********************************************************************/
/*  Vektory pro ulozeni slozek pro transport                        */
/********************************************************************/
#ifdef SKUPINA_TRANSPORT
	EX double   *P_rslo;
	EX double   *P_rslo_new;
 	 #ifdef SKUPINA_DUALPOROSITY
	EX double  *P_rslo_por;
 	 #endif // SKUPINA_DUALPOROSITY
	EX double   *P_sod;
#endif // SKUPINA_TRANSPORT
/********************************************************************/
/*                      Definice datovych struktur                  */
/********************************************************************/
/*   Globalni parametry modelu                                      */
/********************************************************************/
struct S_glp
{
    char        jmeno_MMF[ MAXPATH ];
/*------------------------------------------------------------------*/
    char        jmeno_site[ MAXPATH ];
/*------------------------------------------------------------------*/
    int     nmuzl;                 /* pocet multiuzlu rovinnych     */
    int     nuzl;                  /* pocet uzlu prostorovych       */
    int     nmelm;                 /* pocet multielm rovinnych      */
    int     nmkoef;                /* pocet koeficientu pro MElm    */
    int     nelm;                  /* pocet elementu prostorovych   */
    int     nelm_0;                /* pocet placatych elementu      */
    int     nkoef;                 /* pocet koeficientu pro element */
    int     nsmume;                /* pocet sousedu MU->ME          */
/*------------------------------------------------------------------*/
    int     nstnex;                /* pocet vnejsich sten           */
    int     nstnin;                /* pocet vnitrnich sten          */
/*------------------------------------------------------------------*/
   int      nvrst;                 /* pocet vrstev v siti           */
/*------------------------------------------------------------------*/
   int      nsez;                  /* pocet kroku scenare vypoctu   */
   int      noke;                  /* pocet radku OKE               */
/*------------------------------------------------------------------*/
#ifdef SKUPINY_MATERIALU
   int      nmatr;                 /* pocet skupin materialu        */
#endif // SKUPINY_MATERIALU
/*------------------------------------------------------------------*/
#ifdef SKUPINA_DUALPOROSITY
   BOOL     sw_DP_apply;       /* pocitat s dvoji porozitou     */
   int      ndpor;                 /* pocet skupin DUAL_POROSITY    */
#endif // SKUPINA_DUALPOROSITY
/*------------------------------------------------------------------*/
#ifdef SKUPINA_INTERTRANSFER
   BOOL     sw_CTR_apply;      /* pocitat INTERTRANSFER         */
   int      nctr;                  /* pocet skupin INTERTRANSFER    */
#endif // SKUPINA_INTERTRANSFER
/*------------------------------------------------------------------*/
#ifdef SKUPINA_TRANSPORT
    int     nvod;                  /* pocet typu vod                */
/*------------------------------------------------------------------*/
	BOOL   time_analysis;         /* rozbor casoveho kroku         */
	BOOL    quazi_unsteady;        /* kvazi nestacionarni uloha     */
    BOOL    write_waters;          /* vypsat nalezene vody          */
/*..................................................................*/
    char    jmeno_rini[ MAXPATH ];
    char    jmeno_rpop[ MAXPATH ];
#ifdef SKUPINA_DUALPOROSITY
    char    jmeno_sini[ MAXPATH ];  // mh: zvlast POC pro slepe pory
    char    jmeno_spop[ MAXPATH ];  // mh: zvlast POP pro slepe pory
#endif // SKUPINA_DUALPOROSITY
    int     nrslo;                 /* pocet slozek Roztoku - POP    */
/*..................................................................*/
    int     nrlimit;               /* pocet limitu Roztoku          */
    char    plimit[ NLENSPO +1 ];  /* popis limitni slozky Roztoku  */
    double  *rlimit;               /* hodnoty limitu   Roztoku      */
    int     irlimit;               /* index limitni slozky Roztoku  */
/*..................................................................*/
    double  Eps_Qs;        /* Maximalni chyba bilance vnitrni steny */
    double  Eps_Qe;        /* Maximalni chyba bilance Elementu      */
#endif // SKUPINA_TRANSPORT
/*------------------------------------------------------------------*/
#ifdef SKUPINA_CHEMIE
   BOOL     sw_Chem_WriteSet; /* Vypsat nastaveni chemie  0-NE, 1-ano */
   BOOL     sw_Chem_Consist;  /* Vyzadovat konzistenci    0-ne, 1-ANO */
   char     jmeno_chem_cfg[ MAXPATH ];
   char     jmeno_chem_err[ MAXPATH ];
   TChemie  *Chem_Chem;
   TChemieBloku  **Chem_Blok;
   char     jmeno_hini[ MAXPATH ];
   char     jmeno_hpop[ MAXPATH ];
   int      nhslo;                  /* pocet slozek Horniny - POP   */
   BOOL     sw_kapalna;             /* Pocita se s KAPALNOU fazi    */
   BOOL     sw_pevna;               /* Pocita se s PEVNOU fazi      */
   BOOL     sw_plyn;                /* Pocita se s PLYNNOU fazi     */
//   BOOL     sw_sorpce;              /* Pocita se s SORBOVANOU fazi  */
   BOOL     sw_redox;               /* Pocita se s REDOX fazi       */
   BOOL     sw_radio;               /* Pocita se s RADIO fazi       */
   BOOL     sw_vlastnosti;          /* Pocita se s VLASTNOSTMI      */
#endif // SKUPINA_CHEMIE
/*------------------------------------------------------------------*/
    int     pis_log;
/*------------------------------------------------------------------*/
    char        sit_popis[ 80 ];    /* popis site                   */
    float       sit_alfa;           /* uhel anizotropie             */
   int      sit_nvrst;              /* max. pocet vrstev site       */
};
EX struct S_glp G_glp;
/********************************************************************/
/*   Popis MultiElementu a Elementu                 ( soubor .STE ) */
/********************************************************************/
struct S_melm                     // multielement
{
    int        oznac;             /* oznaceni Melementu  < 0 zrusen  */
    int        ipelm;             /* inter. cislo prvniho elementu   */
    int        npelm;             /* pocet elementu v tomtom Melm    */
    int        muzl[ 3 ];         /* interni cisla prirazenych Muzlu */
    int        ipz0[ 3 ];         /* poradove cislo startovni z-ovky */
    float      k1lm1[ MAXK1LM1 ]; /* koeficienty MElm                */
    int        oblast;            /* cislo oblasti obsahujici Melm   */
    int        odv;                /**/
    int        dov;                /**/
/*------------------------------------------------------------------*/
    double     sumtok;             /* souhrny prutok MELM do okoli  */
};
EX struct S_melm *P_melm;
/********************************************************************/
struct S_elm                       // element
{
    int         imelm;             /* interni cislo Multielementu   */
    int         ivrst;             /* cislo vrtvy                   */
    int         imatr;             /* cislo materialu               */
    		// sem je nactena hodnota z *.stm, je kopirovano do idpor (pro dvoji porozitu, jinak se nepouziva)
    float       koef[ MAXELKOEF ]; /* koeficienty elementu          */
/*------------------------------------------------------------------*/
#ifdef SKUPINA_DUALPOROSITY
    int         idpor;             /* index S_DPOR                  */
#endif // SKUPINA_DUALPOROSITY
/*------------------------------------------------------------------*/
    double      vyska;             /* piezometricka vyska  v T ELM  */
    double      tlak;              /* tlak v T ELM                  */
    double      stntlk[ 5 ];       /* tlak v T steny                */
    double      stntok[ 5 ];       /* pretok stenou                 */
    double      bilance;           /* bilance ELM                   */
    int         stnvod[ 6 ];       /* typ vody prochazejici stenou  */
/*------------------------------------------------------------------*/
    double      objem;             /* objem ELM                     */
    double      porobjm;           /* porovy objem ELM              */
    double      porobjm_por;       /* porovy objem ELM poru         */
/*------------------------------------------------------------------*/
#ifdef SKUPINA_TRANSPORT
    double      *rslo;             /* koncentrace slozek roztoku    */
    double      *rslonew;          /* koncentrace slozek roztoku    */
 #ifdef SKUPINA_DUALPOROSITY
    double      *rslo_por;         /* koncentrace slozek poroveho roztoku */
    double      *rslonew_por;      /* koncentrace slozek poroveho roztoku */
 #endif // SKUPINA_DUALPOROSITY
#endif // SKUPINA_TRANSPORT
/*------------------------------------------------------------------*/
#ifdef SKUPINA_CHEMIE
    TChemieElementuD    *Chem;
    double      *vslo;             /* slozky vlastnost              */
    double      *hslo;             /* koncentrace slozek horniny    */
    double      *pslo;             /* koncentrace slozek plynny     */
 #ifdef SKUPINA_DUALPOROSITY
    TChemieElementuD    *Chem_por;
    double      *vslo_por;    /* slozky porove vlastnost  */
    double      *hslo_por;    /* koncentrace slozek porove horniny   */
    double      *pslo_por;    /* koncentrace slozek porove plynny    */
 #endif // SKUPINA_DUALPOROSITY
#endif // SKUPINA_CHEMIE
/*------------------------------------------------------------------*/
#ifdef SKUPINA_REAKCE
    double      sorpcni_plocha;    /* sorpcni plocha ELM poru       */
    double      *hslo;             /* koncentrace slozek horniny    */
//    double      *sslo;             /* koncentrace slozek sorbovano  */
    double      *pslo;             /* koncentrace slozek plynny     */
 #ifdef SKUPINA_DUALPOROSITY
    double      sorpcni_plocha_por; /* sorpcni plocha ELM porove horniny*/
    double      *hslo_por;    /* koncentrace slozek porove horniny   */
//    double      *sslo_por;    /* koncentrace slozek porove sorbovano */
    double      *pslo_por;    /* koncentrace slozek porove plynny    */
 #endif // SKUPINA_DUALPOROSITY
#endif // SKUPINA_REAKCE
/*------------------------------------------------------------------*/
   double      cas_koef;
   BOOL         cas_typ;
};
EX struct S_elm *P_elm;
/********************************************************************/
/*   Popis MultiUzlu a Uzlu                         ( soubor .STU ) */
/********************************************************************/
struct S_muzl
{
    int     oznac;                    /* oznaceni multiuzlu         */
    int     ismelm;                   /* prvni sousedni multielem   */
    int     nsmelm;                   /* pocet sousednich multielem */
    int     ipuzl;                    /* interni cislo prvniho uzlu */
    int     npuzl;                    /* pocet uzlu v multiuzlu     */
    double  x, y;                     /* souradnice multiuzlu       */
    double  z_povrch;                 /* Z-tova souradnice povrchu  */
    int     odp;                      /**/
    int     dop;                      /**/
};
EX struct S_muzl *P_muzl;
/********************************************************************/
struct S_uzl
{
    int     imuzl;                    /* interni cislo MultiUzlu    */
    double  z0;                       /* vyska daneho uzlu          */
};
EX struct S_uzl *P_uzl;
/********************************************************************/
/*   Popis sousednich MultiElm k danemu MultiUzl                    */
/********************************************************************/
struct S_smume
{
    int muzl;                           /* multiuzel                */
    int melm;                           /* sousedni multielement    */
};
EX struct S_smume *P_smume;
/********************************************************************/
/*   Popis typu cerpanych a vtlacenych vod                          */
/********************************************************************/
struct S_vod
{
    int     typ;                        /* typ vody                 */
    char    nazev[ 128 ];               /* nazev typu vody          */
/*------------------------------------------------------------------*/
	BOOL    sw_OSM;                     /* Zapisovat hmotnosti vy. latek */
    BOOL    sw_OSC;                     /* Zapisovat koncentraci    */
/*------------------------------------------------------------------*/
//  Okamzite koncentrace
    double  *vrslo;       /* koncentrace slozek vtlaceneho ROZTOKU */
    double  *crslo;       /* koncentrace slozek cerpaneho ROZTOKU  */
/*------------------------------------------------------------------*/
//  Absolutni objemy a hmotnosti od casu NULA
    double  Cobjem;                     /* objem cerpane vody             */
    double  *Chmota;                    /* hmotnost cerpanych slozek      */
    double  Vobjem;                     /* objem vtlacene vody            */
    double  *Vhmota;                    /* hmotnost vtlacenych slozek     */
/*------------------------------------------------------------------*/
//  Relativni objemy a hmotnosti od zacatku kroku TRAN
    double  CKobjem;                    /* objem cerpane vody             */
    double  *CKhmota;                   /* hmotnost cerpanych slozek      */
    double  VKobjem;                    /* objem vtlacene vody            */
    double  *VKhmota;                   /* hmotnost vtlacenych slozek     */
/*------------------------------------------------------------------*/
//  Pracovni objemy a hmotnosti
    double  WCobjem;                    /* objem cerpane vody             */
    double  *WChmota;                   /* hmotnost cerpanych slozek      */
    double  WVobjem;                    /* objem vtlacene vody            */
    double  *WVhmota;                   /* hmotnost vtlacenych slozek     */
};
#ifdef SKUPINA_TRANSPORT
EX struct S_vod *P_vod;
#endif // SKUPINA_TRANSPORT
/********************************************************************/
/*   Popis a obsah vrstev                                           */
/********************************************************************/
struct S_pvr
{
    char        popis[ NLENPVR +1 ];    /* popis (oznaceni) vrstvy   */
    int     ipslo;                      /* interni cislo prvni slozky  */
};
EX struct S_pvr *P_pvr;
/********************************************************************/
/*   Popis slozek Roztoku                                           */
/********************************************************************/
struct S_rpo
{
   char        popis[ NLENSPO +1 ];     /* popis slozky Roztoku           */
   char     unit[ 5 ];                  /* jednotka slozky [g/l mg/l]     */
   int      to_chem;                    /* ukazatel do vektoru chemie     */
#ifdef SKUPINA_DUALPOROSITY
	// mh:
   double	difus_koef_DP;				/* relativni koef difuze pro dualni porozitu*/
#endif // SKUPINA_DUALPOROSITY

/*------------------------------------------------------------------*/
};
#ifdef SKUPINA_TRANSPORT
	EX struct S_rpo *P_rpo;
#endif // SKUPINA_TRANSPORT
/********************************************************************/
/*   Popis slozek Horniny                                           */
/********************************************************************/
struct S_hpo
{
   char     popis[ NLENSPO +1 ];        /* popis slozky Horniny           */
   char     unit[ 4 ];                  /* jednotka slozky [hm%]          */
   int      to_chem;                    /* ukazatel do vektoru chemie     */
/*------------------------------------------------------------------*/
};
#ifdef SKUPINA_CHEMIE
	EX struct S_hpo *P_hpo;
#endif // SKUPINA_CHEMIE
/********************************************************************/
/*   Seznam okrajovych podminek                                     */
/********************************************************************/
struct S_sez
{
    char    fname[ MAXPATH ];
    char    txt[ 80 ];                  /* Popis dane OKP                 */
    int     ipoke;                      /* interni cislo OKE              */
    int     npoke;                      /* Pocet OKE                      */
    double  srazky;                     /* Destove srazky   v mm . rok^-1 */
    double  dt;                         /* Casovy krok                    */
    int     NK_Flow;                    /* Pocet casovych kroku FLOW      */
/*------------------------------------------------------------------*/
#ifdef SKUPINA_TRANSPORT
    char    iname[ MAXPATH ];           /* Jmeno INI souboru OKP          */
    int     NK_Tran;                    /* Pocet casovych kroku TRAN      */
#endif // SKUPINA_TRANSPORT
/*------------------------------------------------------------------*/
#ifdef SKUPINA_CHEMIE
    int     NK_Chem;                    /* Pocet casovych kroku CHEM      */
#endif // SKUPINA_CHEMIE
/*------------------------------------------------------------------*/
    int     NK_Result;						// pocet casovych kroku pro vystup do TS3 a BIN
    int     NK_sez;                     /* Pocet casovych kroku F_T_Ch    */
};
#ifdef SKUPINA_SCENAR
	EX struct S_sez *P_sez;
#endif // SKUPINA_SCENAR
/********************************************************************/
/*   Popis okrajovych podminek - novy typ OKE                       */
/********************************************************************/
struct S_oke
{
    int     typ;                        /* typ okrajove podminky          */
    int     ivoda;                      /* typ vody                       */
    float   hodnota_oke;                /* zadana hodnota   OKE           */
    float   koeficient_oke;             /* koeficient OKE                 */
/*------------------------------------------------------------------*/
    int      ie_od;                     /* interni cislo Elm    (dolniho) */
    int      ie_do;                     /* interni cislo Elm    (horniho) */
    int     iis;                        /* interni cislo steny            */
/*------------------------------------------------------------------*/
    int     zapnuto;                    /* priznak zapnuti OKE            */
};
#ifdef SKUPINA_SCENAR
	EX struct S_oke *P_oke;
#endif // SKUPINA_SCENAR
/********************************************************************/
/*   Popis typu materialu pro nenasycene proudeni (z INI souboru)   */
/********************************************************************/
struct S_matr
{
   int      typ;                        /* typ materialu                  */
   int      funkce;                     /* funkce pro dany material       */
   double   S_min;                      /* Minimalni saturace             */
   double   Kr_min;                     /* Minimalni koef. rel. prop.     */
   int      nkoef;                      /* pocet koef dane funkce         */
   double   koef[ MAXMATRKOEF ];        /* koeficienty pro danou funkci   */
};
#ifdef SKUPINY_MATERIALU
	EX struct S_matr *P_matr;
#endif // SKUPINY_MATERIALU
/********************************************************************/
/*   Popis typu materialu pro DUAL POROSITY (z INI souboru)         */
/********************************************************************/
struct S_dpor
{
   int      typ;                        /* typ materialu                */
   int      funkce;                     /* funkce pro dany material     */
   int      nkoef;                      /* pocet koef. dane funkce      */
   double   koef[ MAXDPORKOEF ];        /* koeficienty pro danou funkci */
};
#ifdef SKUPINA_DUALPOROSITY
	EX struct S_dpor *P_dpor;
#endif // SKUPINA_DUALPOROSITY
/********************************************************************/
/*   Popis prechodu pro INTERTRANSFER (z INI souboru)               */
/********************************************************************/
struct S_ctr
{
    int     ivrst;                     /* index sledovane vrstvy    */
   int      ikoef_melm;                /* index prepocivaciho koef. */
   double   UP_hmota;
   double   DOWN_hmota;
};
#ifdef SKUPINA_INTERTRANSFER
	EX struct S_ctr *P_ctr;
#endif // SKUPINA_INTERTRANSFER
//
/********************************************************************/
/*                      Definice funkci                             */
/********************************************************************/
/*  Cteni vstupnich dat                                             */
/********************************************************************/
void ctimmf( void );            /* Nacteni parametru ze souboru   .INI  */
void ctiuzl( void );            /* Nacteni uzlu ze souboru        .UZL  */
void ctielm( void );            /* Nacteni elementu ze souboru    .ELM  */
void ctistu( void );            /* Nacteni uzlu ze souboru        .STU  */
void ctiste( void );            /* Nacteni elementu ze souboru    .STE  */
void ctistm( void );            /* Nacteni koef. Elm ze souboru   .STM  */
void smume( void );             /* Vytvoreni struktury SMUME            */
/*------------------------------------------------------------------*/
void ctihdm( int, int );        /*  Nacteni souboru hydrodynamiky .HDM  */
void nastav_slozky( void );
void uvolni_slozky( void );
void ctipop( int );             /*  Nacteni pocatecnich podminek  .POP  */
/*------------------------------------------------------------------*/
void inicializace_S_vod( void );
void uvolneni_S_vod( void );
/*------------------------------------------------------------------*/
void sumace( double );

//mh:
#ifdef SKUPINA_DUALPOROSITY
	void init_DP_mater( void );  // nastaveni indexu materialu pro DP
#endif // SKUPINA_DUALPOROSITY

/********************************************************************/
/*  Zapis vystupnich dat                                            */
/********************************************************************/
void pis_TS2( double );         /* Zapis souboru roztoku po vrst.  .TS2  */
#ifdef SKUPINA_DUALPOROSITY
	int pisTS34( int, int, double ); /* Zapis vysledku do .TS3 nebo     .TS4  */
 #else //SKUPINA_DUALPOROSITY
	int pisTS34( int, double );     /* Zapis vysledku do .TS3 nebo     .TS4  */
#endif // SKUPINA_DUALPOROSITY
int pisDF0( double, double, int ); /* Zapis vysledku do            .DF0  */
int pisDF1( double, double, int ); /* Zapis vysledku do            .DF1  */
int pisDF2( double, double, int ); /* Zapis vysledku do            .DF2  */
// mh: (doplneno df3 na vystup koncentraci)
int pisDF3( double, double, int ); /* Zapis vysledku do		  .DF3  */
void otevri_bin( int );
void pis_bin( int, float * );
void zavri_bin( void );
int pisPOP(char*) ;  // zapis POP z hodnot na konci vypoctu
/********************************************************************/
/*  Pripravne vypocty                                               */
/********************************************************************/
void *seznam_sten( int );       /*  Vytvari seznam sousednich sten  */
/********************************************************************/
/*  Utility pro pripravne vypocty                                   */
/********************************************************************/
int  uzl_pro_elm( int, int );   /* vraci cislo uzlu pro ELM a iUZL  */
int  iuzl_pro_elm( int, int );  /*  vraci iUZL pro ELM a UZL        */
void buble_sort( int *, int );  /*  Bublinove trideni               */
/********************************************************************/
/*  Vlastni vypocet                                                 */
/********************************************************************/
void vypocet( void );               /* Ridici program vypoctu       */
void sestmat( void );               /* Sestaveni globalni matice    */
void objemy( void );                /* vypocet objemu   elementu    */
void nej_spolky( int, int, int, int *, int * );
int  caskrok( double );             /* vypocet casoveho kroku       */
/********************************************************************/
/*  Souborove UTILITY                                               */
/********************************************************************/
char *fjmeno( char *, char * );         /* Cele jmenu souboru       */
FILE *fotevri( char *, const char * );  /* Otevreni souboru         */
int  aktualnejsi_soubor( char *, char * );
/********************************************************************/
/*  Behove UTILITY                                                  */
/********************************************************************/
void Start_programu( char * );      /* Pocatecni hlavicka programu     */
void Konec_programu( char * );      /* Hlavicka ukonceni behu program  */
void Prerus_program( char *, int ); /* Hlavicka ukonceni behu programu */
/********************************************************************/
/*  Casove UTILITY                                                  */
/********************************************************************/
char *mezi_cas( int );              /* Vraci string - mezicas vypoctu  */
/********************************************************************/
/*  Pametove UTILITY                                                */
/********************************************************************/
void malo_pameti( char *, char *, int );       /* Hlaska malo pameti */
void pole_je_male( int, char *, int, char * ); /* Hlaska male pole  */
/********************************************************************/
/*  Logovaci UTILITY                                                */
/********************************************************************/
void smazlog( void );
void pislog( char *, ... );
void pisscr( char *, ... );
/********************************************************************/
/*  Vystup na obrazovku pro TRANSPORT                               */
/********************************************************************/
void win_tran_START( char * );
void win_tran_STOP( void );
void pis_Tran( char *fmt, ... );
/*------------------------------------------------------------------*/
/*                                          Testy                   */
/*------------------------------------------------------------------*/
void cesta_k_souboru_sestav( char *, char *, char * );
/********************************************************************/
/*  UTILITY pro sledovani verzi                                     */
/********************************************************************/
void main_verze( void );
void ctimmf_verze( struct S_verze * );
void vypocet_verze( struct S_verze * );
void df_file_verze( struct S_verze * );
void pisdf0_verze( struct S_verze * );
void pisdf1_verze( struct S_verze * );
void pisdf2_verze( struct S_verze * );
void s_vody_verze( struct S_verze * );
void caskrok_verze( struct S_verze * );
void win_tran_verze( struct S_verze * );
void slozky_verze( struct S_verze * );
/********************************************************************/
/*  copied from che_head.h to simplify structure of inclusions      */
/********************************************************************/
/* vim:  set ts=3 sw=3 expandtab: */

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

/********************************************************************/
/*  copied from cti_ichnew.c to simplify structure of inclusions    */
/********************************************************************/

float che_poradi(int param1, double param2, double param3);
void ctiich_obecne(void);
void ctiich_latkyvefazi(void);
void ctiich_dalsilatky(void);
void ctiich_reakce(void);

/********************************************************************/
/*  copied from semchem_interface.hh to simplify structure of inclusions    */
/********************************************************************/
void che_outpocp_soubor(FILE *fw);
void che_pocitej_soubor(char *soubor, int *poc_krok);
void che_vypis_soubor(char *soubor);
void che_presun_poc_p_(void);
void che_vypis__soubor(char *soubor);

#endif
