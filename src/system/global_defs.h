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
 * @file    global_defs.h
 * @brief   Global macros to enhance readability and debugging, general constants.
 */

#ifndef GLOBAL_DEFS_H
#define GLOBAL_DEFS_H

//Make sure that headers necessary for following macros are included.
#include <stdint.h>
#include <string.h>
#include <stdlib.h>

#include "system/exc_common.hh"
#include "config.h"
#include "mpi.h"
#include "logger.hh"

/*! @brief Debugging macros.
 *
 *  The macro OLD_ASSERT has to be used for assertion tests. An error occures if
 *  given condition is violated.  Macro accepts additional variables to print.
 *
 *  Example:
 *  @verbatim
 *  OLD_ASSERT( i<size , "Array X overflow: index %d >= alocated size %d.\n",i,size);
 *  @endverbatim
 *
 *  The macro INPUT_CHECK should be used for assertions about user input. So
 *  they produce User Error instead of Program error.
 *
 *  The macro DBGCOUT should be used for debugging messages,
 *  so they can be removed in production version.
 *
 *  OLD_WARN_ASSERT - can be used for consistency tests in debugging version.
 *
 *  @{
 */
 #define INPUT_CHECK(i,...)   do { if (!(i))   xprintf(UsrErr,__VA_ARGS__); } while (0)

/**
 * Actually there are following debugging switches
 * FLOW123D_DEBUG_MESSAGES  - use various debugging messages introduced by DBGCOUT
 * FLOW123D_DEBUG_ASSERTS - use assertion checks introduced by ASSERT
 * FLOW123D_DEBUG_PROFILER - use profiling introduced by START_TIMER, END_TIMER
 *
 * You can turn all off defining: FLOW123D_NODEBUG
 * or turn all on defining: FLOW123D_DEBUG
 *
 * FLOW123D_DEBUG overrides FLOW123D_NODEBUG
 */


#ifdef FLOW123D_NODEBUG

#undef  FLOW123D_DEBUG_MESSAGES
#undef  FLOW123D_DEBUG_ASSERTS
#undef  FLOW123D_DEBUG_PROFILER

#endif


#ifdef FLOW123D_DEBUG

#define  FLOW123D_DEBUG_MESSAGES
#define  FLOW123D_DEBUG_ASSERTS
#define  FLOW123D_DEBUG_PROFILER

#endif


#include "asserts.hh"



#ifdef FLOW123D_DEBUG_ASSERTS

/**
 * Just quick hack to fix some unit tests.
 * TODO:
 * We should make better implementation, rather minimizing
 * usage of macros. And make robust "system" part, that
 * is MPI aware, but not MPI dependent.
 */
//#ifdef FLOW123D_DEBUG_ASSERTS_WITHOUT_MPI
//#define MPI_Comm_rank(A, B)
//#endif // FLOW123D_DEBUG_ASSERTS_WITHOUT_MPI

#define OLD_ASSERT(i,...)   do {\
    if (!(i))  {\
        char msg[1024];\
        sprintf( msg, __VA_ARGS__);\
        int rank=-1;\
        THROW( ExcAssertMsg() << EI_Message(std::string(msg)) << EI_MPI_Rank(rank) );\
    }} while (0)

//#define OLD_WARN_ASSERT(i,...) do { if (!(i))    xprintf(Warn,__VA_ARGS__); } while (0)

#define OLD_ASSERT_EQUAL( a, b)  do {\
    stringstream ss; ss << (a) << " != " << (b); \
    OLD_ASSERT( ((a) == (b)), "Violated assert: %s == %s,\n observed: %s.\n",#a, #b, ss.str().c_str()); \
    } while (0)

#define OLD_ASSERT_LESS( a, b) do {\
    stringstream ss; ss << (a) << " >= " << (b); \
    OLD_ASSERT( ((a) < (b)) , "Violated assert: %s < %s,\n observed: %s.\n",#a,#b, ss.str().c_str()); \
    } while (0)

#define OLD_ASSERT_LE( a, b) do {\
    stringstream ss; ss << (a) << " > " << (b); \
    OLD_ASSERT( ((a) <= (b)) , "Violated assert: %s <= %s,\n observed: %s.\n",#a,#b, ss.str().c_str()); \
    } while (0)

#define OLD_ASSERT_PTR( ptr ) do {\
    OLD_ASSERT( ((ptr) != nullptr) , "Null pointer: %s\n", #ptr ); \
    } while (0)
#else

#define OLD_ASSERT(...)
#define OLD_WARN_ASSERT(...)
#define OLD_ASSERT_EQUAL( a, b)
#define OLD_ASSERT_LESS( a, b)
#define OLD_ASSERT_LE( a, b)
#define OLD_ASSERT_PTR( ptr )

#endif




#ifdef FLOW123D_DEBUG_MESSAGES

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

#define DBGCOUT(...)
#define DBGVAR(var)

#endif




/**
 * These macros are necessary in classes that contain Input::Type::Abstract (PARENT macro)
 * and in classes contain descendant of this Abstract (CHILD macro) if these descendants
 * are initialized through methods of @p Input::Factory class.
 *
 * These macros are necessary for initializing of static variables in classes that contain
 * descendants of parent Abstract.
 */
#define FLOW123D_FORCE_LINK_IN_CHILD(x) int force_link_##x = 0;
#define FLOW123D_FORCE_LINK_IN_PARENT(x) extern int force_link_##x; void func_##x() { force_link_##x = 1; }


///@}

#endif // GLOBAL_DEFS_H
