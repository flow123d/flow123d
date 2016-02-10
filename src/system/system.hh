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
 * @file    system.hh
 * @brief   
 */

#ifndef SYSTEM_H
#define SYSTEM_H


#include <mpi.h>
#include <iostream>

#include "system/global_defs.h"
#include "system/exc_common.hh"


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


// Assuming all compilers supports CXX11 features
#define OPERATOR_NEW_THROW_EXCEPTION



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


char * 	get_log_fname( void );
char * 	get_log_file( void );
void	resume_log_file( void );


#define xprintf(...) _xprintf(__FILE__, __func__, __LINE__, __VA_ARGS__)

int     _xprintf(const char * const xprintf_file, const char * const xprintf_func, const int xprintf_line, MessageType type, const char * const fmt, ... );
void *	xmalloc(size_t size);
void * xrealloc( void * ptr, size_t size );

// TODO: implement as a templated function
#ifndef xfree
    #define xfree(p) \
    do { if (p) { free(p); (p)=NULL; } \
    } while (0) /// test & free memory
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

void *operator new (std::size_t size) OPERATOR_NEW_THROW_EXCEPTION;
void *operator new[] (std::size_t size) OPERATOR_NEW_THROW_EXCEPTION;
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
//! @}

// string operations
char * xstrcpy(const char*);
char * xstrtok(char *s, int position = -1);
char * xstrtok(char*,const char* delim, int position = -1);
int    xchomp( char * s );

/**
 * Wrapper to check return codes of C functions. In particular PETSC calls.
 */
inline void chkerr(unsigned int ierr) {
	do {
		if (ierr != 0) THROW( ExcChkErr() << EI_ErrCode(ierr));
	} while (0);
}

/**
 * Wrapper to check return codes of C functions. In particular PETSC calls.
 * This version do the check just as an debugging assert. So the code is empty
 * in release version.
 */
inline void chkerr_assert(unsigned int ierr) {
#ifdef FLOW123D_DEBUG_ASSERTS
    chkerr(ierr);
#endif
}

#endif
//-----------------------------------------------------------------------------
// vim: set cindent:
