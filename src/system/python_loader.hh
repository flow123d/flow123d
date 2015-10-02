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
 * @file    python_loader.hh
 * @brief   
 */

#ifndef PYTHON_UTILS_HH_
#define PYTHON_UTILS_HH_

#include "global_defs.h"

#ifdef FLOW123D_HAVE_PYTHON

#include "Python.h"
#include <string>


namespace internal {
/**
 * Class to implement initialization and finalization of the python module
 */
class PythonRunning {
public:
    PythonRunning(const std::string &python_home);
    ~PythonRunning();
};
} // close namespace internal



/**
 * Class with static only members, should be used to load and compile Python sources either from file or from string.
 * Implement correct initialization and finalization.
 * TODO:
 * - check validity of results and throw exceptions
 * - better loading of modules from files
 * - explain and resolve RuntimeWarning:
 *   "flow123d_python_loader:1: RuntimeWarning: Parent module 'field_python_script' not found while handling absolute import"
 *   that appears during field_python_test.cpp
 *
 */
class PythonLoader {
public:

	/**
	 * Definition of exception thrown by python compiler or interpreter.
	 */
	TYPEDEF_ERR_INFO(EI_PythonMessage, std::string);
    DECLARE_EXCEPTION(ExcPythonError,
            << "Python Error: " << EI_PythonMessage::val << "\n");


    /**
     * Calls python initialization and guarantee that appropriate finalization will be called.
     * Do nothing if initialization was done.
     *
     * The method with no parameters is called at the beginning of every function of this class, so an explicit call
     * to it has only sense if one would like to provide alternative python home directories.
     * The string has form <prefix>[:<exec_prefix>] where <prefix> is prefix for platform independent libraries
     * (namely sources of python standard libraries) and <exec_prefix> is prefix for platform dependent libraries namely
     * for the python executable.
     */

    static void initialize(const std::string &python_home="");

    /**
     * This function loads a module from the given file.
     * Resulting module has to be deallocated by Py_DECREF() macro.
     */
    static PyObject * load_module_from_file(const std::string& fname);
    /**
     * This function compile code in the given string and creates a module with name @p module_name.
     * Resulting module has to be deallocated by Py_DECREF() macro.
     */
    static PyObject * load_module_from_string(const std::string& module_name, const std::string& source_string);
    /**
     * Method which loads module by given module_name
     * module_name can (and probably will) contain packages path (will contain dots '.' which detonates package)
     *
     * Example:
     *   PyObject * python_module = PythonLoader::load_module_by_name ("profiler.profiler_formatter_module")
     *
     * will import module 'profiler_formatter_module' from package 'profiler'
     */
    static PyObject * load_module_by_name(const std::string& module_name);
    /**
     * Tests whether the error indicator is set, if yes formats and throws exception.
     */
    static void check_error();
    /**
     * Check if python function is callable, if not throws exception.
     */
    static PyObject * get_callable(PyObject *module, const std::string &func_name);
};



#endif // FLOW123D_HAVE_PYTHON

#endif /* PYTHON_UTILS_HH_ */
