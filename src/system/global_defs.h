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
 * @file   global_defs.h
 * @brief  Global macros to enhance readability and debugging, general constants.
 *
 */

#ifndef GLOBAL_DEFS_H
#define GLOBAL_DEFS_H

//Make sure that headers necessary for following macros are included.
#include <stdint.h>
#include <string.h>
#include <stdlib.h>

#include "system/exc_common.hh"


/// @brief Global constants.
/// @{

/// size of input line buffer at various places
#define LINE_SIZE 65536

/// @}


/// @brief Macros to enhance readability
/// @{

#define NDEF  -1    ///< not defined positive integer (obsolete - ints should be initialized by value)
#define NONULL(p)   ((p) != NULL)   /// true for non-null pointer


// set array of pointers of given size to NULL
#define SET_ARRAY_NULL(ptrs,n) do {int i; for(i=0;i<n;i++) *(ptrs+i)=NULL;} while (0)
// set array of ints, float, or doubles to 0
#define SET_ARRAY_ZERO(array,n) memset((array), 0, sizeof(*(array))*(n))

// Useful for printing boolean ints
#define BOOL_2_STR(i) 	( (i) ? "No" : "Yes" )

/*! @brief Debugging macros.
 *
 *  The macro ASSERT has to be used for assertion tests. An error occures if
 *  given condition is violated.  Macro accepts additional variables to print.
 *
 *  Example:
 *  @verbatim
 *  ASSERT( i<size , "Array X overflow: index %d >= alocated size %d.\n",i,size);
 *  @endverbatim
 *
 *  The macro INPUT_CHECK should be used for assertions about user input. So
 *  they produce User Error instead of Program error.
 *
 *  The macro DBGMSG should be used for debugging messages,
 *  so they can be removed in production version.
 *
 *  WARN_ASSERT - can be used for consistency tests in debugging version.
 *
 *  @{
 */
#define INPUT_CHECK(i,...)   do { if (!(i))   xprintf(UsrErr,__VA_ARGS__); } while (0)

/**
 * Actually there are following debugging switches
 * DEBUG_MESSAGES  - use various debugging messages introduced by DBGMSG
 * DEBUG_ASSERTS - use assertion checks introduced by ASSERT
 * DEBUG_PROFILER - use profiling introduced by START_TIMER, END_TIMER
 * DEBUG_FUNCTION_STACK  - use function stack introduced by F_ENTRY
 *
 * You can turn all off defining: Flow123d_NODEBUG
 * or turn all on defining: Flow123d_DEBUG
 *
 * Flow123d_DEBUG overrides Flow123d_NODEBUG
 */

#ifdef Flow123d_NODEBUG

#undef  DEBUG_MESSAGES
#undef  DEBUG_ASSERTS
#undef  DEBUG_PROFILER
#undef  DEBUG_FUNCTION_STACK

#endif


#ifdef Flow123d_DEBUG

#define  DEBUG_MESSAGES
#define  DEBUG_ASSERTS
#define  DEBUG_PROFILER
#define  DEBUG_FUNCTION_STACK

#endif


#ifdef DEBUG_ASSERTS

#define ASSERT(i,...)   do {\
    if (!(i))  {\
        char msg[1024];\
        sprintf( msg, __VA_ARGS__);\
        THROW( ExcAssertMsg() << EI_Message(std::string(msg)) );\
    }} while (0)

#define WARN_ASSERT(i,...) do { if (!(i))    xprintf(Warn,__VA_ARGS__); } while (0)

#else

#define ASSERT(...)
#define WARN_ASSERT(...)

#endif



#ifdef DEBUG_ASSERTS

#define ASSERT_EQUAL( a, b)  do {\
    stringstream ss; ss << (a) << " != " << (b); \
    ASSERT( ((a) == (b)), "Violated assert: %s == %s,\n observed: %s.\n",#a, #b, ss.str().c_str()); \
    } while (0)
#else

#define ASSERT_EQUAL( a, b)

#endif



#ifdef DEBUG_ASSERTS

#define ASSERT_LESS( a, b) do {\
    stringstream ss; ss << (a) << " >= " << (b); \
    ASSERT( ((a) < (b)) , "Violated assert: %s < %s,\n observed: %s.\n",#a,#b, ss.str().c_str()); \
    } while (0)

#else

#define ASSERT_LESS( a, b)

#endif



#ifdef DEBUG_MESSAGES

#define DBGMSG(...) do { xprintf(MsgDbg,__VA_ARGS__); fflush(NULL); } while (0)

/// this is simple macro for debugging output of array of ints
/// Use this with care !!!
#define DBGPRINT_INT(name,size,idx) \
    do {\
        int i__;\
        xprintf(Msg,    "%s (int array size=%d):\n",(name),(size));\
        for(i__=0;i__<(size);i__++) \
            xprintf(Msg,"i: %d int: %d\n",i__,(idx)[i__]); \
    } while (0)




/**
 * Usage:
 * DBGCOUT( << "xy" <<endl );
 */
#define DBGCOUT(...) do { std::cout << "    DBG (" \
									<< __FILE__ << ", " \
									<< __func__ << ", " \
									<< __LINE__ << ")" \
									__VA_ARGS__; } while(0)


/**
 * Usage:
 * DBGVAR( arma::vec( "1 2 3" ) );
 */
#define DBGVAR( var ) DBGCOUT( << #var << " = " << var << endl )

#else

#define DBGMSG(...)
#define DBGPRINT_INT(...)
#define DBGCOUT(...)
#define DBGVAR(var)

#endif

///@}

#endif // GLOBAL_DEFS_H
