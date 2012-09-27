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
 * @brief ???
 *
 */

#include <string>
#include "input/input_type.hh"

using namespace std;

#ifndef MAIN_H
#define MAIN_H



// problem types
/*
typedef enum {
CONVERT_TO_OUTPUT=0,
STEADY_SATURATED=1,
UNSTEADY_SATURATED=2,
PROBLEM_DENSITY=3,
UNSTEADY_SATURATED_LMH=4} ProblemType;
*/


class Application {
public:
    Application(const int argc, char ** argv);
    static Input::Type::Record &get_input_type();
    void parse_cmd_line(const int argc, char ** argv);
    void free_and_exit();
private:
    string main_input_dir_;
    string main_input_filename_;
    string log_filename_;
    int passed_argc_;
    char ** passed_argv_;
};




#endif

//-----------------------------------------------------------------------------
// vim: set cindent:

