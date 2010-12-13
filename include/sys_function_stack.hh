/*!
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
#include <sstream>
#include <string>
#include <list>


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

#ifdef NODEBUG

#define F_STACK_SHOW(stream)
#define F_ENTRY
#define F_ENTRY_P(param_string)

#else

/*
 *  __SOMETHING__ stands for name of local instance of the Trace class. That means that
 *  F_ENTRY macro can be used only once in each function (otherwise it is redeclaration error).
 */
#define F_STACK_SHOW(stream) flow::Trace::stack_print(stream)
#define F_ENTRY              flow::Trace __SOMETHING__( __FILE__, __func__, __LINE__ , "")
#define F_ENTRY_P(param_string) flow::Trace __SOMETHING__( __FILE__, __func__, __LINE__ , param_string )

#endif

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
     */
    class Trace
    {
    public:
        /**
         * Constructor, takes source filename and func_name of a function, and line numebr line_no
         * of the call of #F_ENTRY. Forms a message string and put it onto stack.
         */
        explicit Trace( const char * filename, const char * func_name, const int line_no , std::string param_string)
        {
            std::ostringstream oss;
            oss << "Tracepoint: " << filename << ", " << func_name << "(), line " << line_no << ", params: " << param_string;
            program_stack.push_back( oss.str() );
        }

        /**
         * Destructor, automatically called on the function exit. Pop out the stack message of that function.
         */
        ~Trace()
        {
            program_stack.pop_back();
        }

        /**
         * Prints the stack into a given file.
         */
        void static stack_print( FILE * fw )
        {
            fprintf( fw, "Stack trace, depth: %u\n", (unsigned int) program_stack.size());
            for( std::list<std::string>::reverse_iterator psi=program_stack.rbegin(); psi!=program_stack.rend(); psi++ )
                fprintf( fw, " %s\n", (*psi).c_str() );
        }

        /**
         * As previous, but output into C++ output stream.
         */
        void static stack_print( std::ostream * os )
        {
            *os << "Stack trace, depth: " << program_stack.size() << std::endl;
            for( std::list<std::string>::iterator psi=program_stack.begin(); psi!=program_stack.end(); psi++ )
                *os << *psi << std::endl;
        }
    private:
        /**
         * Static stack of trace messages.
         */
        static std::list<std::string> program_stack;
    };
}

#endif /* SYS_FUNCTION_STACK_HH_ */
