/*
 * exceptions.cc
 *
 *  Created on: May 9, 2012
 *      Author: jb
 */

#include "exceptions.hh"
#include <iostream>
#include <boost/exception/diagnostic_information.hpp>

/**
 * Error header.
 *
 */
void ExceptionBase::print_exc_data (std::ostream &out) const
{
}

void ExceptionBase::print_info(std::ostream &out) const {
    out << "No error message." << std::endl;
}

const char * ExceptionBase::what() const throw () {
    // Be sure that this function do not throw.
    try {
        std::ostringstream converter;

        converter << "--------------------------------------------------------" << std::endl;


        print_info(converter);
  
        converter << "\n** Diagnosting info **\n" ;
        converter << boost::diagnostic_information_what( *this );
        
        //    print_stack_trace (converter);
        converter << "--------------------------------------------------------" << std::endl;

        // creates local string on stack, we hope that there is a space
        // freed by stack unrolling
        return converter.str().c_str();

    } catch (std::exception &exc) {
        std::cerr << "*** Exception encountered in exception handling routines ***" << std::endl << "*** Message is " << std::endl
                << exc.what() << std::endl << "*** Aborting! ***" << std::endl;

        std::abort();
    } catch (...) {
        std::cerr << "*** Exception encountered in exception handling routines ***" << std::endl << "*** Aborting! ***"
                << std::endl;

        std::abort();
    }
    return 0;
}
