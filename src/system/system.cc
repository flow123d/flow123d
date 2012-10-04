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
 * @ingroup system
 * @brief  Various system-wide functions
 *
 */

#include <cstring>
#include <cstdarg>
#include <ctime>
#include <cstdlib>
#include <sys/stat.h>
#include <cerrno>
#include <sstream>

#include <fstream>
#include <string>

#include "global_defs.h"
#include "system/system.hh"
//#include "io/read_ini.h"
#include "system/xio.h"
#include "system/file_path.hh"

#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/trim.hpp>




SystemInfo sys_info;

static bool petsc_initialized = false;

/*!
 * @brief Read system parameters, open log.
 */
void system_init( MPI_Comm comm,const  string &log_filename )
{
    int ierr;

    //for(int i=0;i<argc;i++) xprintf(Msg,"%s,",argv[i]);
     petsc_initialized = true;
    sys_info.comm=comm; 


    xio_init(); //Initialize XIO library

    // TODO : otevrit docasne log file jeste pred ctenim vstupu (kvuli zachyceni chyb), po nacteni dokoncit
    // inicializaci systemu

    ierr=MPI_Comm_rank(comm, &(sys_info.my_proc));
    ierr+=MPI_Comm_size(comm, &(sys_info.n_proc));
    ASSERT( ierr == MPI_SUCCESS,"MPI not initialized.\n");

    DBGMSG("MPI size: %d rank: %d\n",sys_info.n_proc,sys_info.my_proc);

    // determine logfile name or switch it off
    stringstream log_name;

    if (log_filename != "") {
        if ( log_filename == "\n" ) {
           // -l option without given name -> turn logging off
           sys_info.log=NULL;
    } else {
      // given log name
           log_name << log_filename <<  "." << sys_info.my_proc << ".log";
           sys_info.log_fname = FilePath(log_name.str(), FilePath::output_file );
           sys_info.log=xfopen(sys_info.log_fname.c_str(),"wt");

        }
    } else {
        // use default name
        log_name << "flow123."<< sys_info.my_proc << ".log";
        sys_info.log_fname = FilePath(log_name.str(), FilePath::output_file );
        sys_info.log=xfopen(sys_info.log_fname.c_str(),"wt");

    }

    sys_info.verbosity=0;
    sys_info.pause_after_run=0;
}

void system_set_from_options()
{
    //sys_info.verbosity = OptGetInt( "Run", "Screen_verbosity", "0" );

}


/// @brief INTERNAL DEFINITIONS FOR XPRINTF
/// @{

struct MsgFmt {
	int  num;           ///< format number
	bool log;		    ///< log the message - YES/NO
	bool mpi;           ///< treat as global message (invoke MPI_Barrier() when printing)
	int screen;	        ///< print to stdout,stderr,NULL
    bool stop;          ///< terminate the program
	const char * head;	///< message formating string
};

#define SCR_NONE	0
#define SCR_STDOUT	1
#define SCR_STDERR	2

/// configuration table for individual message types defined in system.h
/// Msg type    Log    mpi      screen      Stop    message header
#define	NUM_OF_FMTS		8
static struct MsgFmt msg_fmt[] = {
	{Msg, 		true,  false,   SCR_STDOUT,	false,	NULL},
	{MsgDbg,    true,  false,   SCR_STDOUT, false,  "DBG (%s, %s(), %d):"},
	{MsgLog,	true,  false,   SCR_NONE,	false,	NULL},
	{MsgVerb,	false, false,   SCR_STDOUT,	false,	NULL},
	{Warn,		true,  false,   SCR_STDERR,	false,	"Warning (%s, %s(), %d):\n"},
	{UsrErr,	true,  false,   SCR_STDERR,	true,	"User Error (%s, %s(), %d):\n"},
	{Err,		true,  false,   SCR_STDERR,	true,	"Error (%s, %s(), %d):\n"},
	{PrgErr,	true,  false,   SCR_STDERR, true,	"Internal Error (%s, %s(), %d):\n"}
};

/// @}

/*!
 * @brief Multi-purpose printing routine: messages, warnings, errors
 * @param[in] xprintf_file   current file
 * @param[in] xprintf_func   current function
 * @param[in] xprintf_line   current line number
 * @param[in] type           message type
 * @param[in] fmt            message format
 * @return      Same as printf, what internal printing routine returns.
 */
int _xprintf(const char * const xprintf_file, const char * const xprintf_func, const int xprintf_line, MessageType type, const char * const fmt, ... )
{
	struct MsgFmt	mf;
	int rc;
	FILE *screen=NULL;

	static int mpi_msg_id = 0;
	int ierr;

	if ((type == MsgVerb) && (sys_info.verbosity <= 0))
		return 0;

	if ((type < 0) || (type >= NUM_OF_FMTS))
		type = Msg;
	mf = msg_fmt[type];

	// determine output stream
	switch (mf.screen) {
		case SCR_STDOUT : screen=stdout; break;
		case SCR_STDERR : screen=stderr; break;
		case SCR_NONE : screen=NULL; break;
		default: screen=NULL;
	}


#ifdef MPI_PRINTALL
	//print msg header with mpi_id to distinguish the message origin
    if (screen)
        fprintf(screen,"[mpi_id=%d] ", sys_info.my_proc );
    if (mf.log && sys_info.log)
        fprintf(sys_info.log,"[mpi_id=%d] ", sys_info.my_proc );
#else
    //if not PRINTALL, allow console output only for MPI master, no need to print mpi_id
	if ( (screen) && ( sys_info.my_proc != 0 ) )
        screen = NULL;
#endif

    //generate barrier and unique ID for MPI messages
    if (mf.mpi) {
        ierr = MPI_Barrier(sys_info.comm);
        if (ierr != MPI_SUCCESS ) {
            printf("MPI_Barrier() error in xprintf()\n"); //can not call xprintf() when xprintf() is failing
            exit(EXIT_FAILURE);
        }

        // print global msg_id
        if (screen)
            fprintf(screen,"[msg_id=%d] ", mpi_msg_id );
        if (mf.log && sys_info.log)
            fprintf(sys_info.log,"[msg_id=%d] ", mpi_msg_id );
        mpi_msg_id++;
    }

	// print head
	if (mf.head) {
		if (screen) fprintf(screen,mf.head,xprintf_file,xprintf_func,xprintf_line);
		if (mf.log && sys_info.log)
			fprintf(sys_info.log,mf.head,xprintf_file,xprintf_func,xprintf_line);
	}
	// print message
	{
		va_list argptr;

		if (mf.log && sys_info.log)
		{
			va_start( argptr, fmt );
			rc=vfprintf(sys_info.log,fmt,argptr); //rc=char written, <0 if err
			va_end( argptr );
            // flush every message (maybe there is a problem in cygwin without that)
		    fflush(sys_info.log);
		}

		if (screen)
		{
			va_start( argptr, fmt );
			rc=vfprintf(screen,fmt,argptr);
			va_end( argptr );
		    // flush every message (maybe there is a problem in cygwin without that)
			fflush(screen);
		}
	}


	if (mf.stop) {
	    // explicit flush of all streams
		fflush(NULL);
        xterminate(true);
	}
	return rc;
}

/*!
 * @brief     Terminates the program.
 *
 * @param[in] on_error Set true if the function is called as a result of an error. This produce backtrace.
 *
 * TODO: should be destructor of a main program object with pointer to the main application object, through deleting
 * application object it should delete all created objects, namely free all memory and close all files.
 * 
 * More over the application should be derived from ApplicationBase which collects functionality of system.cc
 * In descendants  of ApplicationBase we can call PetscFinalize, and keep ApplicationBase free of petsc.
 *
 */

#ifdef HAVE_PETSC

#include <petsc.h>
#endif

int xterminate( bool on_error )
{
    //close the Profiler
    Profiler::uninitialize();


    if (on_error) { F_STACK_SHOW( stderr ); }

	//TODO: Free memory, close files

#ifdef HAVE_PETSC	
	if ( petsc_initialized )
	{
           PetscErrorCode ierr=0;
 
	   ierr = PetscFinalize(); CHKERRQ(ierr);
           
           on_error = (ierr != 0);
	   petsc_initialized = false;
	}
#endif
    if (sys_info.log) xfclose( sys_info.log );

    if (sys_info.pause_after_run) {
        printf("\nPress <ENTER> for closing the window\n");
        getchar();
    }

    //fflush and fclose all files (including stdout, stderr, stdio)
    //this function is GNU extension
    fcloseall();

    //select proper Return Code
    if ( on_error ) //error in program or during PetscFinalize()
        exit( EXIT_FAILURE );
    else
        exit( EXIT_SUCCESS );
}

/*!
 * @brief Memory allocation with checking.
 *
 * Allocates memory block with checking of correct size and successful allocation.
 *
 * @param[in] size  New size for the memory block, in bytes.
 * @return          same as ISO C realloc()
 */
void *xmalloc( size_t size )
{
	void *rc;


	if (size == 0 ) size++;
	//ASSERT( size != 0 ,"Bad requested size (%u bytes) for memory allocation\n", size );
	rc = malloc( size );
	if ( rc == NULL ) xprintf(Err ,"Not enough memory for allocating %u bytes\n", size );

	return(rc);
}

/*!
 * @brief Reallocation of memory block with checking.
 *
 * Reallocates memory block with checking of successful reallocation.
 * The size of the memory block pointed to by the ptr parameter is changed to the size bytes,
 * expanding or reducing the amount of memory available in the block.
 *
 * The function may move the memory block to a new location, in which case the new location
 * is returned. The content of the memory block is preserved up to the lesser of
 * the new and old sizes, even if the block is moved. If the new size is larger, the value
 * of the newly allocated portion is indeterminate.
 *
 * In case that ptr is NULL, the function behaves exactly as malloc,
 * assigning a new block of size bytes and returning a pointer to the beginning of it.
 *
 * In case that the size is 0, the memory previously allocated in ptr is deallocated
 * as if a call to free was made, and a NULL pointer is returned.
 *
 * @param[in] ptr   Pointer to a memory block previously allocated
 * @param[in] size  New size for the memory block, in bytes.
 * @return          same as ISO C realloc()
 */
void * xrealloc( void * ptr, size_t size )
{
    void * rc;

    F_ENTRY;

    rc = realloc( ptr, size );
    if ( rc == NULL ) xprintf(Err ,"Not enough memory for allocating %u bytes\n", size );

    return(rc);
}

/*
void *operator new (std::size_t size, const my_new_t &) throw() {
    return xmalloc(size);
}

void *operator new[] (std::size_t size, const my_new_t &) throw() {
    return xmalloc(size);
}

void operator delete( void *p,  const my_new_t &) throw ()
{
    xfree(p);
}

void operator delete[]( void *p,  const my_new_t &) throw ()
{
    xfree(p);
}
*/

void *operator new (std::size_t size) throw(std::bad_alloc) {
    return xmalloc(size);
}

void *operator new[] (std::size_t size) throw(std::bad_alloc) {
    return xmalloc(size);
}

void operator delete( void *p) throw()
{
    xfree(p);
}

void operator delete[]( void *p) throw()
{
    xfree(p);
}



/*!
 * @brief SYSTEM with err handling
 */
int xsystem( const char *cmd )
{
	int rc;

	F_ENTRY;

	rc = system( cmd );
	INPUT_CHECK(!( rc != 0 ),"Error executing external command: %s\n", cmd );
	return(rc);
}

/*!
 * @brief MAKE BRAND NEW COPY OF STRING
 */
char *xstrcpy( const char *src )
{
	char *rc;
	size_t length;

	F_ENTRY;

	ASSERT(!( src == NULL ),"NULL pointer as argument of function xstrcpy()\n");
	length = strlen( src ) + 1;
	rc = (char*) xmalloc(length * sizeof(char));
	strcpy( rc, src );
	return(rc);
}

/*!
 * @brief      STRTOK WITH ERROR HANDLING and whitespace delimiters
 * @param[in]  s        strtok string pointer
 * @param[in]  position requested position of the token
 * @return              strtok return
 */
char *xstrtok(char *s, int position)
{
    char *rc;
    const char * const whitespace_delim=" \t\r\n";

    F_ENTRY;

    rc = xstrtok( s, whitespace_delim, position);
    return(rc);
}

/*!
 * @brief      STRTOK WITH ERROR HANDLING and user specified delimiters
 * @param[in]  s1       strtok string pointer
 * @param[in]  delim    delimiters
 * @param[in]  position requested position of the token
 * @return              strtok return
 *
 * Function behaves like original strtok
 */
char *xstrtok( char *s1, const char *delim, int position )
{
	char *rc;
	static char * full_string = NULL;
	static int token_count;

	F_ENTRY;

	ASSERT(!( delim == NULL ),"NULL pointer as delimiter in xstrtok()\n");

	if ( s1 )
	{
	    if ( !full_string )
	    {
	        full_string = (char *)xmalloc( LINE_SIZE );
	        full_string[0] = 0x0;
	    }

	    strncpy( full_string, s1, LINE_SIZE );
	    token_count = 0;
    }

    INPUT_CHECK( token_count == position || position < 0, "Requested position %d dosn't match token position %d", position, token_count);
	rc = strtok( s1, delim );
	token_count++;

	INPUT_CHECK(!( rc == NULL ),"Missing token no. %d: original string '%s' with delimiters '%s'\n", token_count, full_string, delim );

	return(rc);
}

/*!
 * @brief           Delete trailing whitespace characters (space,tab,CR,NL).
 * @param[in,out] s string to change
 * @return          number of deleted characters
 */
int xchomp( char * s )
{
    int no_erased = 0;
    char * p;

    F_ENTRY;

    ASSERT(NONULL(s), "Can not chomp NULL string.");

    if ( *s ) //string not empty
    {
        p = s;
        while (*p)
            p++; //find end of string
        p--; //set p to the last character
        while ((p >= s) && ((*p == ' ') || (*p == '\t') || (*p == '\r') || (*p == '\n')))
        {
            *p = 0x0;
            no_erased++;
            p--;
        }
    }
    return(no_erased);
}


/*!
 * @brief MKDIR WITH ERROR HANDLING
 */
int xmkdir( const char *s )
{
    int rc;

    F_ENTRY;

    ASSERT(!( s == NULL ),"NULL pointer as argument of function xmkdir()\n");
    rc = mkdir(s, S_IRWXU); // create dir with rwx perm. for user
    if (errno == EEXIST)
        rc = 0;
    INPUT_CHECK(!( rc != 0 ),"Cannot make directory %s\n", s );

    return(rc);
}

/*!
 * @brief RMDIR with err handling
 */
int xrmdir( const char *s )
{
    int rc;

    F_ENTRY;

    ASSERT(!( s == NULL ),"NULL pointer as argument of function xrmdir()\n");
    rc = rmdir( s );
    INPUT_CHECK(!( rc != 0 ),"Cannot delete directory %s\n", s );
    return(rc);
}

/*!
 * @brief CHDIR WITH ERROR HANDLING
 */
int xchdir( const char *s )
{
    int rc;

    F_ENTRY;

    ASSERT(!( s == NULL ),"NULL pointer as argument of function xchdir()\n");
    rc = chdir( s );
    INPUT_CHECK(!( rc != 0 ),"Cannot change directory to %s\n", s );
    return(rc);
}

/*!
 * @brief DELETE a FILE with error handling
 */
int xremove( const char *fname )
{
    int rc;

    F_ENTRY;

    ASSERT(!( fname == NULL ),"NULL pointer as argument of function xremove()\n");
    if( access( fname , F_OK ) == 0 )
    {
        rc = remove( fname );
        INPUT_CHECK(!( rc != 0 ),"Cannot remove file %s\n", fname );
    }
    else
        xprintf( Warn, "File '%s' does not exist, can not remove. Ignoring.\n", fname );

    return(rc);
}

/*!
 * @brief GET CURRENT WORKING DIRECTORY with error handling
 */
char *xgetcwd( void )
{
    char tmp[PATH_MAX];
    char * rc;

    F_ENTRY;

    rc = getcwd( tmp, PATH_MAX );
    ASSERT(NONULL(rc),"Cannot get name of current working directory\n");

    return(xstrcpy( tmp ));
}

/*!
 *  @brief Skip to given section in a given file.
 *  @param[in,out]  in          Handle of the file that we search.
 *  @param[in]      section     Section name to find.
 *  @return                     true - if we have found the section, false otherwise
 */
bool skip_to( FILE *const in, const char *section )
{
    char line[ LINE_SIZE ];
    char string[ LINE_SIZE ];

    F_ENTRY;

    ASSERT( NONULL( in ), "Null input file handle.\n");
    ASSERT( NONULL( section ), "NULL section.\n");

    while( xfgets( line, LINE_SIZE - 2, in ) != NULL ) {
        sscanf( line, "%s", string ); // strip spaces
        if( strcmpi( string, section ) == 0 )
        {
            return( true );
        }
    }

    return(false);
}



/*!
 *  @brief Skip to the first line match  @p pattern up to surrounding spaces and case.
 *  @param[in,out]  in          Input stream to search.
 *  @param[in]      pattern     String to look for.
 *  @return                     true - if we have found the section, false otherwise
 */
bool skip_to( istream &in, const string &pattern )
{
    char line[ LINE_SIZE ];
    char string[ LINE_SIZE ];

    F_ENTRY;

    for(std::string line; ! in.eof() ; std::getline(in, line) ) {
        boost::trim(line);
        if ( boost::iequals( line, pattern ) ) return true;
    }

    return false;
}


