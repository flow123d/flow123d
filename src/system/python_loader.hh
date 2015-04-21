/*
 * python_utils.hh
 *
 *  Created on: Aug 31, 2012
 *      Author: jb
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
