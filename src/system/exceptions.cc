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
 * @file    exceptions.cc
 * @brief   
 */

#include "system.hh"
#include "exceptions.hh"
#include <iostream>
#include <cstring>
#include <boost/exception/diagnostic_information.hpp>


ExceptionBase::ExceptionBase() {
	this->frames_to_cut_ = { "boost", "exception_detail", "throw_exception"};
}



ExceptionBase::ExceptionBase(const ExceptionBase &other)
: stack_trace_(other.stack_trace_), frames_to_cut_(other.frames_to_cut_) {}




ExceptionBase::~ExceptionBase() throw () {}




void ExceptionBase::print_stacktrace(std::ostream &out) const {
	stack_trace_.print(out, frames_to_cut_);
}


const char * ExceptionBase::what() const throw () {
    // have preallocated some space for error message we want to return
    // Is there any difference, if we move this into ExceptionBase ??
    static std::string message(1024,' ');

    // Be sure that this function do not throw.
    try {
        std::ostringstream converter;

        converter << std::endl << std::endl;
        converter << "--------------------------------------------------------" << std::endl;
        converter << this->what_type_msg();
        print_info(converter);

        converter << "\n** Diagnosting info **\n" ;
        converter << boost::diagnostic_information_what( *this );
        print_stacktrace(converter);
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


std::string ExceptionBase::what_type_msg() const {
	return "Program Error: ";
}



