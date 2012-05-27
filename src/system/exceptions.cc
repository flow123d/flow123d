/*
 * exceptions.cc
 *
 *  Created on: May 9, 2012
 *      Author: jb
 */

#include "system.hh"
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

void ExceptionBase::print_info(std::ostringstream &out) const {
    DBGMSG("\n");
    out << "No error message." << std::endl;
}

const char * ExceptionBase::what() const throw () {
    // have preallocated some space for error message we want to return
    // Is there any difference, if we move this into ExceptionBase ??
    static std::string message(1024,' ');

    // Be sure that this function do not throw.
    try {
        std::ostringstream converter;

        converter << std::endl;
        converter << "--------------------------------------------------------" << std::endl;
        converter << "Program Error: ";
        print_info(converter);
        converter << "\n** Diagnosting info **\n" ;
        converter << boost::diagnostic_information_what( *this );
        //    print_stack_trace (converter);
        converter << "--------------------------------------------------------" << std::endl;

        message = converter.str();
        return message.c_str();

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



const char * InputException::what() const throw () {
    // have preallocated some space for error message we want to return
    // Is there any difference, if we move this into ExceptionBase ??
    static std::string message(1024,' ');


    // Be sure that this function do not throw.
    try {
        std::ostringstream converter;

        converter << "--------------------------------------------------------" << std::endl;
        converter << "User Error: ";
        print_info(converter);
        converter << "--------------------------------------------------------" << std::endl;

        message = converter.str();
        return message.c_str();

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
