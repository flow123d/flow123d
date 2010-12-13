/*!
 *
 * $Id$
 * $Revision$
 * $LastChangedBy$
 * $LastChangedDate$
 *
 * @file
 * @brief  Implementation of function call stack.
 *
 */

#include <string>
#include <list>

#include "sys_function_stack.hh"

namespace flow
{
///definition of global variable for stacktrace storage
    std::list<std::string> Trace::program_stack;
}
