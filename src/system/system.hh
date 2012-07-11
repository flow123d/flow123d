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
#include "io/read_ini.h"
#include "system/math_fce.h"
#include "system/sys_function_stack.hh"
#include "system/sys_profiler.hh"
#include "system/xio.h"


// for a linux system we assume glibc library
// with support of ISOC99 functions
//#define _ISOC99_SOURCE
#ifndef _BSD_SOURCE
   #define _BSD_SOURCE
#endif

#define strcmpi strcasecmp
#define PATH_SEP "/"

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

} SystemInfo;

extern SystemInfo sys_info;

void 	system_init( int argc, char ** argv, const string& log_filename);


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
bool    skip_to( FILE *const in, const char *section );
//! @}

// string operations
char * xstrcpy(const char*);
char * xstrtok(char *s, int position = -1);
char * xstrtok(char*,const char* delim, int position = -1);
int    xchomp( char * s );

/**
 * @brief Singleton class used to change given 'variables' in the file path.
 *
 * IONameHandler stores items formed by the combination of a key value and a mapped value.
 * The pairs are used to substitution in paths of the input and output files.
 *
 * If program is running with "-i" switch the "${INPUT}" variable in the given filename is substituted by the value immediately following this switch.
 * If program is running with "-o" switch all output file paths are prefixed by the value immediately following this switch (instead of root directory flow.ini file)
 *
 * @par Example usage:
 * @code
 * IONameHandler::get_instance()->get_input_file_name("relative/path/with/${VAR}/flow.ini");
 * @endcode
 *
 * @par Another example usage:
 * @code
 * IONameHandler &io_name_handler = *(IONameHandler::get_instance());
 * io_name_handler.get_input_file_name("relative/path/with/${INPUT}/var/mesh.msh");
 * @endcode
 */
class IONameHandler {
public:
    /*!
	 * @brief Returns instance of IONameHandler.
	 *
	 * Class IONameHandler is created as a singleton (Singleton design pattern). Static method get_instance() returns pointer to IONameHandler object.
	 *
	 * @par Example usage:
	 * @code
	 * IONameHandler& io_name_handler = *(IONameHandler::get_instance());
	 * @endcode
	 */
	static IONameHandler* get_instance();
	/*!
	 * @brief Returns absolute path to given input file.
	 *
	 * @param[in] file_name Filename relatively to root directory.
	 * @return absolute path to input file
	 */
	string get_input_file_name(string file_name);
	/*!
	 * @brief Returns absolute path to given output file.
	 *
	 * @param[in] 	file_name Filename relatively to output directory.
	 * @return 		absolute path to output file
	 */
	string get_output_file_name(string file_name);
	/*!
	* @brief Returns value of output_dir variable.
	*/
	string get_output_dir();
	/*!
	 * @brief Add new item to place holder.
	 *
	 * IONameHandler placeholder is extended by adding a single new item. The item can be used in the name of the input or output file name.
	 *
	 * @par Example usage:
	 * @code
	 * IONameHandler::get_instance()->add_placeholder_item("${SUBST_VAL}", "path/value");
	 * @endcode
	 *
	 * @param[in] 	key Key of new item.
	 * @param[in] 	val Value of new item.
	 * @return 		always true
	 */
	bool add_placeholder_item(string key,string value);
	/*
	 * @brief Removes the key (and its corresponding value) from place holder.
	 *
	 * @param[in] key The key that needs to be removed.
	 * @return The value to which the key had been mapped in place holder, or empty string if the key did not have a mapping.
	 */
//  string remove_placeholder_item(string key);
private:
  static IONameHandler* instance;
  string root_dir;
  string output_dir;
  std::map<string,string> placeholder;

  /*!
   * @brief Private constructor prevents instantiation from other classes
   */
  IONameHandler() {};
  /*!
   * @brief Private copy constructor
   */
  IONameHandler(IONameHandler const&) {};
  /*!
   * @brief Private assignment operator - can never be called
   */
  IONameHandler& operator=(IONameHandler const&) { return (*this);};
  /*!
   * @brief Initialization of root directory - where is located flow.ini file.
   */
  void initialize_root_dir();
  /*!
   * @brief Initialization of output directory - the value of "-o" command line variable if is set, otherwise where is located flow.ini file.
   */
  void initialize_output_dir();
  /*!
   * @brief initialization of placeholder standard pairs - the value of "-i" command line variable for "${INPUT}" if is set.
   */
  void initialize_placeholder();
  /*!
   * @brief Returns value of root_dir variable.
   */
  string get_root_dir();
  /*!
   * @brief Substitute all variable in the given file path.
   */
  string substitute_value(string file);
};

#endif
//-----------------------------------------------------------------------------
// vim: set cindent:
