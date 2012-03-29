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
#include <string>
#include <vector>


/**
 * @def F_ENTRY
 *
 * Unless NODEBUG is defined, it creates instance of the Trace class, providing compile time information
 * about the place of usage.
 *
 * @def F_ENTRY_P( param_string )
 *
 * Allows include any C++ string into function call report.
 * TODO: Because macros can not be overloaded, this has to have different name then the previous one.
 * Moreover this is not type save. Consider use inline functions. But this leads to change from F_ENTRY to F_ENTRY().
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
#define F_STACK_SHOW(stream) flow::Trace::stack_print(stream)
#define F_ENTRY              flow::Trace __SOMETHING__( __FILE__, __func__, __LINE__ , "")
#define F_ENTRY_P(param_string) flow::Trace __SOMETHING__( __FILE__, __func__, __LINE__ , param_string )


#else

#define F_STACK_SHOW(stream)
#define F_ENTRY
#define F_ENTRY_P(param_string)

#endif

#define MAX_STACK_DEPTH 127

namespace flow
{
    /**
     *  @brief This class provides function call stack.
     *
     *  Then the backtrace of the called functions can be reported in the case of an error.
     *  Operations with the stack should be only through #F_ENTRY and #F_STACK_SHOW(stream).
     *  To include a function in the call stack just use macro #F_ENTRY at the very beginning of the function.
     *  #F_ENTRY macro is nonempty only for the debugging version. For the release version, i.e. when NODEBUG
     *  is defined, #F_ENTRY is empty and produce no run-time overhead.
     *
     *  Stack depth is limited to MAX_STACK_DEPTH entries. Additional entries will cause improper stack trace output,
     *  but NOT program failure.
     *
     */
    class Trace
    {
    public:
        /**
         * Constructor, takes source filename and func_name of a function, and line number line_no
         * of the call of #F_ENTRY. Forms a message string and put it onto stack.
         */
        explicit Trace( const char * filename, const char * func_name, const int line_no , std::string param_string)
        {
//BUFSIZE of 16 should be enough for integer to string conversion
            #define BUFSIZE 16
            static char charbuf[BUFSIZE] = {0,};

            if ( stack_depth >= MAX_STACK_DEPTH )
                return;

            program_stack[stack_depth] = "Tracepoint: ";
            program_stack[stack_depth].append(filename);
            program_stack[stack_depth].append(", ");
            program_stack[stack_depth].append(func_name);
            program_stack[stack_depth].append("(), line ");

            snprintf( charbuf, BUFSIZE-1, "%d", line_no );
            program_stack[stack_depth].append( charbuf );

            program_stack[stack_depth].append(", params: ");
            program_stack[stack_depth].append(param_string);

            stack_depth++;
#undef BUFSIZE
        }

        /**
         * Destructor, automatically called on the function exit. Pop out the stack message of that function.
         */
        ~Trace()
        {
            if ( stack_depth > 0 )
                stack_depth--;
        }

        /**
         * Prints the stack into a given file.
         */
        static void  stack_print( FILE * fw )
        {
            fprintf( fw, "Stack trace, depth: %d\n", stack_depth );
            for( int i = stack_depth-1; i >= 0 ; --i )
                fprintf( fw, " %s\n", program_stack[i].c_str() );
        }

        /**
         * As previous, but output into C++ output stream.
         */
        static void  stack_print( std::ostream * os )
        {
            *os << "Stack trace, depth: " << stack_depth << std::endl;
            for( int i = stack_depth-1; i >= 0 ; --i )
                *os << program_stack[i] << std::endl;
        }
    private:
        /**
         * Static stack of trace messages.
         */
        static int stack_depth;
        static std::vector<std::string> program_stack;
    };
}

#endif /* SYS_FUNCTION_STACK_HH_ */
