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
 * @brief Stack trace, using macro F_ENTRY to put a trace point.
 *
 */

#ifndef SYS_FUNCTION_STACK_HH_
#define SYS_FUNCTION_STACK_HH_

#include <cstdio>
#include <iostream>
#include <vector>


/**
 * @def F_ENTRY
 *
 * Unless NODEBUG is defined, it creates instance of the Trace class, providing compile time information
 * about the place of usage.
 *
 * @def F_STACK_SHOW
 *
 * Just shortcut to both static output methods of the Trace class.
 */

#ifdef DEBUG_FUNCTION_STACK

/*
 *  __SOMETHING__ stands for name of local instance of the Trace class. That means that
 *  F_ENTRY macro can be used only once in each function (otherwise it is redeclaration error).
 */
//macros to convert number to quoted string using preprocessor
#define TO_XSTR(s) TO_STR(s)
#define TO_STR(s) #s

// quoted strings in source are fused into one by preprocessor
//__func__ can not be converted to quoted string and fused with preprocessor - it is created later (in compiler stage)
#define F_ENTRY              flow::Trace __SOMETHING__( "Tracepoint: " __FILE__ ", ", __func__,  "(), line " TO_XSTR(__LINE__) "." )

#define F_STACK_SHOW(stream) flow::Trace::stack_print(stream)

#else

#define F_STACK_SHOW(stream)
#define F_ENTRY

#endif

namespace flow
{
    struct Trace_helper {
        const char * str_file;
        const char * str_func;
        const char * str_line;
    };
    /**
     *  @brief This class provides function call stack.
     *
     *  Then the backtrace of the called functions can be reported in the case of an error.
     *  Operations with the stack should be only through #F_ENTRY and #F_STACK_SHOW(stream).
     *  To include a function in the call stack just use macro #F_ENTRY at the very beginning of the function.
     *  #F_ENTRY macro is nonempty only for the debugging version. For the release version, i.e. when NODEBUG
     *  is defined, #F_ENTRY is empty and produce no run-time overhead.
     *
     */
    class Trace
    {
    public:
        /**
         * Constructor, takes source filename and func_name of a function, and line number line_no
         * of the call of #F_ENTRY. Forms a message string and put it onto stack.
         */
        explicit Trace( const char * const trace_str, const char * const trace_str2, const char * const trace_str3 )
        {
            if (stack_depth == -1) {
                stack_depth=0;
                program_stack=new std::vector<Trace_helper>(16);
            }
            (*program_stack)[stack_depth].str_file = trace_str;
            (*program_stack)[stack_depth].str_func = trace_str2;
            (*program_stack)[stack_depth].str_line = trace_str3;
            stack_depth++;
        }

        /**
         * Destructor, automatically called on the function exit. Pop out the stack message of that function.
         */
        ~Trace()
        {
                stack_depth--;
        }

        /**
         * Prints the stack into a given file.
         */
        static void  stack_print( FILE * fw )
        {
            fprintf( fw, "Stack trace, depth: %d\n", stack_depth );
            for( int i = stack_depth-1; i >= 0 ; --i )
                fprintf( fw, " %s%s%s\n", (*program_stack)[i].str_file, (*program_stack)[i].str_func, (*program_stack)[i].str_line );
        }

        /**
         * As previous, but output into C++ output stream.
         */
        static void  stack_print( std::ostream * os )
        {
            *os << "Stack trace, depth: " << stack_depth << std::endl;
            for( int i = stack_depth-1; i >= 0 ; --i )
                *os << (*program_stack)[i].str_file << (*program_stack)[i].str_func << (*program_stack)[i].str_line << std::endl;
        }
    private:
        /**
         * Static stack of trace messages.
         */
        static int stack_depth;
        static std::vector<Trace_helper> *program_stack;
    };
}

#endif /* SYS_FUNCTION_STACK_HH_ */
