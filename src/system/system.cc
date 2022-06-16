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
 * @file    system.cc
 * @ingroup system
 * @brief   Various system-wide functions
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
#include "mpi.h"

#include "system/system.hh"
#include "system/file_path.hh"
#include "system/sys_profiler.hh"

#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/format.hpp>



SystemInfo sys_info;


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




/*!
 * @brief SYSTEM with err handling
 */
/*
int xsystem( const char *cmd )
{
	int rc;

	rc = system( cmd );
	INPUT_CHECK(!( rc != 0 ),"Error executing external command: %s\n", cmd );
	return(rc);
}*/

/*!
 * @brief MAKE BRAND NEW COPY OF STRING
 */
/*
char *xstrcpy( const char *src )
{
	char *rc;
	size_t length;

	ASSERT_PERMANENT_PTR(src).error("NULL pointer as argument of function xstrcpy()\n");
	length = strlen( src ) + 1;
	rc = (char*) xmalloc(length * sizeof(char));
	strcpy( rc, src );
	return(rc);
}*/

/*!
 * @brief      STRTOK WITH ERROR HANDLING and whitespace delimiters
 * @param[in]  s        strtok string pointer
 * @param[in]  position requested position of the token
 * @return              strtok return
 */
/*
char *xstrtok(char *s, int position)
{
    char *rc;
    const char * const whitespace_delim=" \t\r\n";

    rc = xstrtok( s, whitespace_delim, position);
    return(rc);
}
*/

/*!
 * @brief      STRTOK WITH ERROR HANDLING and user specified delimiters
 * @param[in]  s1       strtok string pointer
 * @param[in]  delim    delimiters
 * @param[in]  position requested position of the token
 * @return              strtok return
 *
 * Function behaves like original strtok
 */
/*
char *xstrtok( char *s1, const char *delim, int position )
{
	char *rc;
	static char * full_string = NULL;
	static int token_count;

	ASSERT_PERMANENT_PTR(delim).error("NULL pointer as delimiter in xstrtok()\n");

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
*/

/*!
 * @brief           Delete trailing whitespace characters (space,tab,CR,NL).
 * @param[in,out] s string to change
 * @return          number of deleted characters
 */
/*
int xchomp( char * s )
{
    int no_erased = 0;
    char * p;

    ASSERT_PERMANENT(s).error("Can not chomp NULL string.");

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
}*/


/*!
 * @brief MKDIR WITH ERROR HANDLING
 */
/*
int xmkdir( const char *s )
{
    int rc;

    ASSERT_PERMANENT_PTR(s).error("NULL pointer as argument of function xmkdir()\n");
    rc = mkdir(s, S_IRWXU); // create dir with rwx perm. for user
    if (errno == EEXIST)
        rc = 0;
    INPUT_CHECK(!( rc != 0 ),"Cannot make directory %s\n", s );

    return(rc);
}*/

/*!
 * @brief RMDIR with err handling
 */
/*
int xrmdir( const char *s )
{
    int rc;

    ASSERT_PERMANENT_PTR(s).error("NULL pointer as argument of function xrmdir()\n");
    rc = rmdir( s );
    INPUT_CHECK(!( rc != 0 ),"Cannot delete directory %s\n", s );
    return(rc);
}*/

/*!
 * @brief CHDIR WITH ERROR HANDLING
 */
/*
int xchdir( const char *s )
{
    int rc;

    ASSERT_PERMANENT_PTR(s).error("NULL pointer as argument of function xchdir()\n");
    rc = chdir( s );
    INPUT_CHECK(!( rc != 0 ),"Cannot change directory to %s\n", s );
    return(rc);
}*/

/*!
 * @brief DELETE a FILE with error handling
 */

/*
int xremove( const char *fname )
{
    int rc;

    ASSERT_PERMANENT_PTR(fname).error("NULL pointer as argument of function xremove()\n");
    if( access( fname , F_OK ) == 0 )
    {
        rc = remove( fname );
        INPUT_CHECK(!( rc != 0 ),"Cannot remove file %s\n", fname );
    }
    else
    	WarningOut() << "File '" << fname << "' does not exist, can not remove. Ignoring." << std::endl;

    return(rc);
}*/

/*!
 * @brief GET CURRENT WORKING DIRECTORY with error handling
 */
/*
char *xgetcwd( void )
{
    char tmp[PATH_MAX];
    char * rc;

    rc = getcwd( tmp, PATH_MAX );
    ASSERT_PERMANENT_PTR(rc(.error("Cannot get name of current working directory\n");

    return(xstrcpy( tmp ));
}*/

