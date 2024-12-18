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

/**
 * Actually there are following debugging switches
 * FLOW123D_DEBUG_MESSAGES  - use various debugging messages introduced by DBGCOUT
 * FLOW123D_DEBUG_ASSERTS - use assertion checks introduced by ASSERT_PERMANENT
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
#define _TOKENPASTE(x, y) func_ ## x ## y     // helper macro
#define _TOKENPASTE2(x, y) _TOKENPASTE(x, y)  // helper macro
#define FLOW123D_FORCE_LINK_IN_PARENT(x) extern int force_link_##x; void _TOKENPASTE2(x, __LINE__)(void) { force_link_##x = 1; }


///@}

#endif // GLOBAL_DEFS_H
