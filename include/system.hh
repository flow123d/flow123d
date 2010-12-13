#ifndef SYSTEM_H
#define SYSTEM_H

#include <mpi.h>


#include "global_defs.h"
#include "read_ini.h"
#include "math_fce.h"
#include "sys_function_stack.hh"

// for a linux system we assume glibc library
// with support of ISOC99 functions
//#define _ISOC99_SOURCE
#ifndef _BSD_SOURCE
   #define _BSD_SOURCE
#endif
#include <stdio.h>

#define strcmpi strcasecmp
#define PATH_SEP "/"

using namespace std;


// **************************************************************
/*!  @brief  Identifiers for various output messages.
 */
typedef enum {
    Msg = 0, MsgDbg, MsgLog, MsgVerb, Warn, UsrErr, Err, PrgErr
} MessageType;

// **************************************************************
/*!  @brief  Timing structure for simple timing messages in log file.
 */
#define TIMING_TAG_SIZE 100

typedef struct Timing {
	char tag[TIMING_TAG_SIZE];			// identification signature used in timing log messages
	const char *msg_info;			// "global"/"local"
	double start, last;		// time of call timing_create, last time of call timing_meantime
	MPI_Comm comm;			// PETSC_COMM_WORLD or PETSC_COMM_SELF
	int logging_proc; 		// which proc writes message
} Timing;

// **************************************************************
/*!  @brief  System structure for various global variables.
 */
typedef struct SystemInfo {
    int  verbosity;             // system printf verbosity
    int pause_after_run;        // to keep terminal open on Windows
    char * log_fname;           // name of the master log file
    FILE *log;                  // log file handle
    Timing *timing;             // whole program timing

    int n_proc;                 // number of processors
    int my_proc;                // self processor number

} SystemInfo;

extern SystemInfo sys_info;

void 	system_init( int &argc, char ** &argv);
void    system_set_from_options();
char * 	get_log_fname( void );
char * 	get_log_file( void );
void	resume_log_file( void );

Timing * timing_create(const char * tag,MPI_Comm comm);
void     timing_meantime(Timing * t);
void     timing_destroy(Timing *t);
void     timing_reuse(Timing *t, const char * tag);

#define xprintf(...) _xprintf(__FILE__, __func__, __LINE__, __VA_ARGS__)

int     _xprintf(const char * const xprintf_file, const char * const xprintf_func, const int xprintf_line, MessageType type, const char * const fmt, ... );
int    	xterminate( bool on_error );
void *	xmalloc(size_t size);
void * xrealloc( void * ptr, size_t size );

// TODO: implement as a templated function
#ifndef xfree
    #define xfree(p) do { if (p) { free(p); (p)=NULL; } else {DBGMSG("Free NULL pointer? (in %s, %s(), line %d)\n", __FILE__, __func__, __LINE__);} } while (0) /// test & free memory
#endif


/**
 * @brief Replacement of new/delete operator in the spirit of xmalloc.
 *
 * Up to my knowledge overloading of original new/delete is the only clean.
 * Possibly disadvantage is that all 'new' calls in system and other templates
 * become also overloaded.
 *
 */
// @{

void *operator new (std::size_t size) throw(std::bad_alloc);
void *operator new[] (std::size_t size) throw(std::bad_alloc);
void operator delete( void *p) throw();
void operator delete[]( void *p) throw();
// @}

int     xsystem(const char*);

//! @brief Operations on files and directories
/// @{
int     xmkdir( const char *s );  ///< Create directory (GLIBC function, original in <sys/stat.h>)
int     xchdir( const char *s );  ///< Change directory (GLIBC function, original in <unistd.h>)
int     xremove( const char *s ); ///< Remove file or directory (function)
char *  xgetcwd( void );          ///< Get current working directory (GLIBC function, original in <unistd.h>)
int     xrename( const char * oldname, const char * newname ); ///< Rename file (function)
//#define tmpfile xtmpfile  NOT_USED    ///< Open a temporary file (function)
//! @}

//! @brief auxiliary
/// @{
bool    skip_to( FILE *const in, const char *section );
//! @}

// string operations
char * xstrcpy(const char*);
char * xstrtok(char *s, int position = -1);
char * xstrtok(char*,const char* delim, int position = -1);
int    xchomp( char * s );

#endif
//-----------------------------------------------------------------------------
// vim: set cindent:
