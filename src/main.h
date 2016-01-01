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
 * @file    main.h
 * @brief   
 */

#include <string>
#include "input/input_type_forward.hh"
#include "input/accessors.hh"
#include "input/type_output.hh"
#include "system/application_base.hh"

using namespace std;

#ifndef MAIN_H
#define MAIN_H



class Application : public ApplicationBase {
public:
    /// Root of the Input::Type tree. Description of whole input structure.
    static Input::Type::Record & get_input_type();
    
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

    /**
     * Split path to directory (part up to last DIR_SEPARATOR excluded) and filename.
     * Directory is set to "." if no DIR_SEPARATOR is found.
     */
    void split_path(const string& path, string& directory, string& file_name);

private:

    /// Get version of program and other base data from rev_num.h and store them to map
    Input::Type::RevNumData get_rev_num_data();

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

