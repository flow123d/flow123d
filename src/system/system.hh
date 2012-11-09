/*!
 *
 * Copyright (C) 2007 Technical University of Liberec.  All rights reserved.
 *
 * Please make a following refer to Flow123d on your project site if you use the program for any purpose,
 * especially for academic research:
 * Flow123d, Research Centre: Advanced Remedial Technologies, Technical University of Liberec, Czech Republic
 *
 * This program is free software; you can redistribute it and/or modify it under the terms
 * of the GNU General Public License version 3 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 021110-1307, USA.
 *
 *
 * $Id$
 * $Revision$
 * $LastChangedBy$
 * $LastChangedDate$
 *
 * @file
 * @brief ???
 *
 */

#ifndef SYSTEM_H
#define SYSTEM_H

//#include <mpi>
#include <cstdio>
#include <string>


#include "global_defs.h"
#include "system/math_fce.h"
#include "system/sys_function_stack.hh"
#include "mpi.h"


// for a linux system we assume glibc library
// with support of ISOC99 functions
//#define _ISOC99_SOURCE
#ifndef _BSD_SOURCE
   #define _BSD_SOURCE
#endif

#ifdef WINDOWS
  #include <direct.h>
  #define GetCurrentDir _getcwd
  #define PATH_MAX MAX_PATH 
#else
  #include <unistd.h>
  #define GetCurrentDir getcwd
  #include <limits.h>   //PATH_MAX
#endif

#define strcmpi strcasecmp
#define DIR_DELIMITER '/'

using namespace std;


// **************************************************************
/*!  @brief  Identifiers for various output messages.
 */
typedef enum MessageType {
    Msg = 0, MsgDbg, MsgLog, MsgVerb, Warn, UsrErr, Err, PrgErr
} MessageType;


// **************************************************************
/*!  @brief  System structure for various global variables.
 */
typedef struct SystemInfo {
    int  verbosity;             // system printf verbosity
    int pause_after_run;        // to keep terminal open on Windows
    string log_fname;           // name of the master log file
    FILE *log;                  // log file handle

    int n_proc;                 // number of processors
    int my_proc;                // self processor number
    
    MPI_Comm comm;

} SystemInfo;

extern SystemInfo sys_info;

void system_init( MPI_Comm comm,const  string &log_fname );


void    system_set_from_options();
char * 	get_log_fname( void );
char * 	get_log_file( void );
void	resume_log_file( void );


#define xprintf(...) _xprintf(__FILE__, __func__, __LINE__, __VA_ARGS__)

int     _xprintf(const char * const xprintf_file, const char * const xprintf_func, const int xprintf_line, MessageType type, const char * const fmt, ... );
int    	xterminate( bool on_error );
void *	xmalloc(size_t size);
void * xrealloc( void * ptr, size_t size );

// TODO: implement as a templated function
#ifndef xfree
    #define xfree(p) \
    do { if (p) { free(p); (p)=NULL; } \
         else {DBGMSG("Free NULL pointer? (in %s, %s(), line %d)\n", __FILE__, __func__, __LINE__); \
              } \
    } while (0) /// test & free memory
#endif
//        F_STACK_SHOW( stdout ); \

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
bool skip_to( FILE *const in, const char *section );
bool skip_to( istream &in, const string &pattern );
//! @}

// string operations
char * xstrcpy(const char*);
char * xstrtok(char *s, int position = -1);
char * xstrtok(char*,const char* delim, int position = -1);
int    xchomp( char * s );
#endif
//-----------------------------------------------------------------------------
// vim: set cindent:
