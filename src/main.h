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
#include "input/accessors.hh"
#include "system/application_base.hh"

using namespace std;

#ifndef MAIN_H
#define MAIN_H



class Application : public ApplicationBase {
public:
    /// Root of the Input::Type tree. Description of whole input structure.
    static Input::Type::Record input_type;
    
    /// Application constructor. 
    Application(int argc, char ** argv);
    
    /**
     * Displays program version and build info.
     * Pass version information to Profiler.
     * 
     * TODO: Split these two functionalities.
     */ 
    void display_version();
    
    /**
     * Read main input file
     * 
     * Returns accessor to the root Record.
     */ 
    Input::Record read_input();
    
    /// Destructor
    virtual ~Application();

protected:

    /**
     * Run application.
     *
     * Read input and solve problem.
     */
    virtual void run();

    /**
     * Check pause_after_run flag defined in input file.
     */
    virtual void after_run();

    /**
     * Parse command line parameters.
     * @param[in] argc       command line argument count
     * @param[in] argv       command line arguments
     */
    virtual void parse_cmd_line(const int argc, char ** argv);

private:

    /// directory of main input file (used to resolve relative paths of other input files)
    string main_input_dir_;
    /// filename of main input file
    string main_input_filename_;

    int passed_argc_;
    char ** passed_argv_;
    
    /// Description of possible command line arguments.
    string program_arguments_desc_;

    /// If true, we do output of profiling information.
    bool use_profiler;

    /// root input record
    Input::Record root_record;
};




#endif

//-----------------------------------------------------------------------------
// vim: set cindent:

