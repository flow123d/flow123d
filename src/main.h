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
#include "system/exceptions.hh"
class HC_ExplicitSequential;


using namespace std;

#ifndef MAIN_H
#define MAIN_H



class Application : public ApplicationBase {
public:
    TYPEDEF_ERR_INFO( EI_InputVersionStr, string);
    DECLARE_EXCEPTION( ExcVersionFormat,
            << "Wrong format of the version specification: "
            << EI_InputVersionStr::qval);
    DECLARE_EXCEPTION( ExcUnknownProblem, << "Problem type not implemented.\n" );


    /// Root of the Input::Type tree. Description of whole input structure.
    static Input::Type::Record & get_input_type();
    
    /// Application constructor. 
    Application(const std::string &python_path);
    
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
    
    /**
     * Run application.
     *
     * Read input and solve problem.
     */
    void run() override;

    /**
     * Terminate all MPI processes if exception is thrown.
     */
    void terminate();

    /// Destructor
    virtual ~Application();

protected:

    /**
     * Check pause_after_run flag defined in input file.
     */
    void after_run();


    /**
     * Parse command line parameters.
     * @param[in] argc       command line argument count
     * @param[in] argv       command line arguments
     */
    virtual void parse_cmd_line(const int argc, char ** argv);

private:

    /// Get version of program and other base data from rev_num.h and store them to map
    Input::Type::RevNumData get_rev_num_data();

    /// Main Flow123d problem
    HC_ExplicitSequential *problem_;

    /// filename of main input file
    string main_input_filename_;

    //int passed_argc_;
    //char ** passed_argv_;
    
    /// Description of possible command line arguments.
    string program_arguments_desc_;

    /// If true, we do output of profiling information.
    bool use_profiler;

    /// location of the profiler report file
    string profiler_path;

    /// If true, preserves output of balance in YAML format.
    bool yaml_balance_output_;

    /// root input record
    Input::Record root_record;
};




#endif

//-----------------------------------------------------------------------------
// vim: set cindent:

