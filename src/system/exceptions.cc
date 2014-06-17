/*
 * exceptions.cc
 *
 *  Created on: May 9, 2012
 *      Author: jb
 */

#include "system.hh"
#include "exceptions.hh"
#include <iostream>
#include <cstring>
#include <boost/exception/diagnostic_information.hpp>

#ifdef HAVE_EXEC_INFO
#include <execinfo.h>
#endif

#ifdef HAVE_DEMAGLER
#include <cxxabi.h>
#endif

ExceptionBase::ExceptionBase()
: stacktrace(NULL),n_stacktrace_frames(0)
{
    fill_stacktrace();
}



ExceptionBase::ExceptionBase(const ExceptionBase &other)
: stacktrace(NULL), n_stacktrace_frames( other.n_stacktrace_frames)
{
    int size=0;
    for(int i=0; i<n_stacktrace_frames; i++) size+=strlen(other.stacktrace[i])+1;
    stacktrace = (char **) malloc(size*sizeof(char) +( n_stacktrace_frames+1)*sizeof(stacktrace) );
    char *ptr_dest = (char *)(stacktrace + n_stacktrace_frames);
    for(int i=0; i<n_stacktrace_frames; i++) {
        stacktrace[i]=ptr_dest;
        char * ptr_src = other.stacktrace[i];
        while(*ptr_src!=0) *(ptr_dest++) = *(ptr_src++);
        *(ptr_dest++)=0;
    }
}




ExceptionBase::~ExceptionBase() throw () {
    //DBGMSG("st: %p %p %d\n", this, stacktrace, n_stacktrace_frames);
    if (stacktrace) {
        free(stacktrace);
        stacktrace=NULL;
        n_stacktrace_frames=0;
        //DBGMSG("st: %p %p %d\n", this, stacktrace, n_stacktrace_frames);
    }
}


void ExceptionBase::fill_stacktrace()
{
#ifdef HAVE_EXEC_INFO
    if (! stacktrace && ! n_stacktrace_frames) {
        void * array[25];
        n_stacktrace_frames = backtrace(array, 25);
        stacktrace = backtrace_symbols(array, n_stacktrace_frames);
        //DBGMSG("st: %p %p %d\n", this, stacktrace, n_stacktrace_frames);
    }
#endif
}








void ExceptionBase::print_stacktrace(std::ostream &out) const {
    out << endl;
    out << "** Stacktrace **" << endl;

    int i_frame;
    for(i_frame=0; i_frame < n_stacktrace_frames; i_frame++) {
        string frame(stacktrace[i_frame]);
        if (   frame.find("boost") != string::npos
            && frame.find("exception_detail") != string::npos
            && frame.find("throw_exception") != string::npos
            ) break;
    }
    i_frame++;

    unsigned int out_i_frame=0;
    for(;i_frame< n_stacktrace_frames; i_frame++, out_i_frame++) {
        string frame(stacktrace[i_frame]);
        unsigned int start_pos = frame.find("(")+1,
                     end_pos = frame.find("+");
        string magled_fname = frame.substr( start_pos, end_pos-start_pos );

        int status=-1;
        char *demagled_f_name = {0};

#ifdef HAVE_DEMAGLER
        demagled_f_name = abi::__cxa_demangle(magled_fname.c_str(), 0, 0, &status);
#endif
        if (status == 0) {
            out << setw(3) << out_i_frame << "  " << demagled_f_name << endl;
            free(demagled_f_name);
        } else {
            out << setw(3) << out_i_frame << "  " << magled_fname << endl;
        }

        if (magled_fname == "main") break;
    }
    out << endl;

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
        converter << "Program Error: ";
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



