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
 * @file    stack_trace.cc
 * @brief
 */

#include "system/stack_trace.hh"
#include "config.h"
#include <boost/exception/diagnostic_information.hpp>

#ifdef FLOW123D_HAVE_EXEC_INFO
#include <execinfo.h>
#endif


#ifdef FLOW123D_HAVE_DEMAGLER
#include <cxxabi.h>
#endif


StackTrace::StackTrace()
: frames_(NULL),n_frames_(0)
{
#ifdef FLOW123D_HAVE_EXEC_INFO
    if (! frames_ && ! n_frames_) {
        void * array[25];
        n_frames_ = backtrace(array, 25);
        frames_ = backtrace_symbols(array, n_frames_);
    }
#endif
}


StackTrace::StackTrace(const StackTrace &other)
: frames_(NULL), n_frames_( other.n_frames_)
{
    int size=0;
    for(int i=0; i<n_frames_; i++) size+=strlen(other.frames_[i])+1;
    frames_ = (char **) malloc(size*sizeof(char) +( n_frames_+1)*sizeof(frames_) );
    char *ptr_dest = (char *)(frames_ + n_frames_);
    for(int i=0; i<n_frames_; i++) {
    	frames_[i]=ptr_dest;
        char * ptr_src = other.frames_[i];
        while(*ptr_src!=0) *(ptr_dest++) = *(ptr_src++);
        *(ptr_dest++)=0;
    }
}


StackTrace::~StackTrace()
{
    if (frames_) {
        free(frames_);
        frames_=NULL;
        n_frames_=0;
    }
}


void StackTrace::print(std::ostream &out, std::vector<std::string> frames_to_cut) const
{
	using namespace std;

    out << endl;
    out << "** Stacktrace **" << endl;

    int i_frame;
    for(i_frame=0; i_frame < n_frames_; i_frame++) {
        string frame(frames_[i_frame]);
        bool is_to_cut = false; // check if frame is intended to cut
        for (auto to_cut : frames_to_cut) {
            if ( frame.find(to_cut) != string::npos ) is_to_cut = true;
        }
        if (is_to_cut) break;
        /*if (   frame.find("boost") != string::npos
            && frame.find("exception_detail") != string::npos
            && frame.find("throw_exception") != string::npos
            ) break;*/
    }
    i_frame++;

    unsigned int out_i_frame=0;
    for(;i_frame< n_frames_; i_frame++, out_i_frame++) {
        string frame(frames_[i_frame]);

        /*bool is_to_cut = false; // check if frame is intended to cut
        for (auto to_cut : frames_to_cut) {
        	if ( frame.find(to_cut) != string::npos ) is_to_cut = true;
        }
        if (is_to_cut) continue;*/

        unsigned int start_pos = frame.find("(")+1,
                     end_pos = frame.find("+");
        string magled_fname = frame.substr( start_pos, end_pos-start_pos );

        int status=-1;
        char *demagled_f_name = {0};

#ifdef FLOW123D_HAVE_DEMAGLER
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
