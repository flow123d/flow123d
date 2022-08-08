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
 * @file    python_loader.cc
 * @brief   
 */

#include "global_defs.h"

#ifdef FLOW123D_HAVE_PYTHON

#include "Python.h"
#include "system/python_loader.hh"

//#include "system/system.hh"
#include "system/file_path.hh"
#include <string>
//#include <iostream>
#include <fstream>
//#include <sstream>
#include <boost/algorithm/string.hpp>
#include <include/pybind11/pybind11.h>
#include <include/pybind11/embed.h> // everything needed for embedding

using namespace std;
namespace py = pybind11;

string python_sys_path = R"CODE(

def get_paths():
    import sys
    import os
    # print('\n'.join(sys.path))
    return os.pathsep.join(sys.path)

)CODE";

// default value
string PythonLoader::sys_path = "";


void PythonLoader::initialize(const std::string &python_home)
{
    static internal::PythonRunning _running(python_home);
}


py::module_ PythonLoader::load_module_from_file(const std::string& fname) {
    initialize();

    py::module_ module = py::module_::import(fname.c_str());
    return module;
}



py::module_ PythonLoader::load_module_from_string(const std::string& module_name, const std::string& func_name, const std::string& source_string) {
    initialize();

    py::dict globals = py::globals();
    py::exec(source_string.c_str(), globals, globals);

    py::module_ user_module = py::module::create_extension_module(module_name.c_str(), "", new PyModuleDef());
    user_module.add_object(func_name.c_str(), globals[func_name.c_str()]);
    return user_module;
}

py::module_ PythonLoader::load_module_by_name(const std::string& module_name) {
    initialize();

    // import module by dot separated path and its name
    py::module_ module_object = py::module_::import (module_name.c_str());
    PythonLoader::check_error();

    return module_object;
}


void PythonLoader::check_error() {
    if (PyErr_Occurred()) {
        PyObject *ptype, *pvalue, *ptraceback;
        PyErr_Fetch(&ptype, &pvalue, &ptraceback);
        string error_description = string(PyUnicode_AsUTF8(PyObject_Str(pvalue)));
        
        // clear error indicator
        PyErr_Clear();
    
        /* See if we can get a full traceback */
        PyObject *traceback_module = PyImport_ImportModule("traceback");
    
        string str_traceback;
        if (traceback_module && ptraceback) {
            PyObject *format_tb_method = PyObject_GetAttrString(traceback_module, "format_tb");
            if (format_tb_method) {
                PyObject * traceback_lines = PyObject_CallFunctionObjArgs(format_tb_method, ptraceback, NULL);
                if (traceback_lines) {
                    PyObject * join_str = PyUnicode_FromString("");
                    PyObject * join = PyObject_GetAttrString(join_str, "join");
                    PyObject * message = PyObject_CallFunctionObjArgs(join, traceback_lines, NULL);
                    if (message) {
                        str_traceback = string(PyUnicode_AsUTF8(PyObject_Str(message)));
                    }
                    Py_DECREF(traceback_lines);
                    Py_DECREF(join_str);
                    Py_DECREF(join);
                    Py_DECREF(message);
                }
            }
        }
        
        // get value of python's "sys.path"
        string python_path = PythonLoader::sys_path;
        replace(python_path.begin(), python_path.end(), ':', '\n');
        
        // construct error message
        string py_message =
                   "\nType: " + string(PyUnicode_AsUTF8(PyObject_Str(ptype))) + "\n"
                 + "Message: " + string(PyUnicode_AsUTF8(PyObject_Str(pvalue))) + "\n"
                 + "Traceback: \n" + str_traceback + "\n"
                 + "Paths: " + "\n" + python_path + "\n";
    
        THROW(ExcPythonError() << EI_PythonMessage( py_message ));
    }
}



//PyObject * PythonLoader::get_callable(PyObject *module, const std::string &func_name) {
//    char func_char[func_name.size()+2];
//    strcpy(func_char, func_name.c_str());
//    PyObject * func = PyObject_GetAttrString(module, func_char );
//    PythonLoader::check_error();
//
//    if (! PyCallable_Check(func)) {
//    	stringstream ss;
//    	ss << "Field '" << func_name << "' from the python module: " << PyModule_GetName(module) << " is not callable." << endl;
//    	THROW(ExcPythonError() << EI_PythonMessage( ss.str() ));
//    }
//
//    return func;
//}



wstring to_py_string(const string &str) {
    wchar_t wbuff[ str.size() ];
    size_t wstr_size = mbstowcs( wbuff, str.c_str(), str.size() );
    return wstring( wbuff, wstr_size );
}

string from_py_string(const wstring &wstr) {
    char buff[ wstr.size() ];
    size_t str_size = wcstombs( buff, wstr.c_str(), wstr.size() );
    return string( buff, str_size );
}

#define STR_EXPAND(tok) #tok
#define STR(tok) string(STR_EXPAND(tok))

namespace internal {

PythonRunning::PythonRunning(const std::string& program_name)
{
#ifdef FLOW123D_PYTHON_PREFIX
        static wstring _python_program_name = to_py_string(program_name);
        Py_SetProgramName( &(_python_program_name[0]) );
        wstring full_program_name = Py_GetProgramFullPath();
        // cout << "full program name: " << from_py_string(full_program_name) << std::endl;

        // try to find string "flow123d" from right side of program_name
        // if such a string is not present, we are most likely unit-testing
        // in that case, full_flow_prefix is current dir '.'
        size_t pos = full_program_name.rfind( to_py_string("flow123d") );
        wstring full_flow_prefix = ".";
        if (pos != wstring::npos) {
            full_flow_prefix = full_program_name.substr(0,pos-string("/bin/").size() );
        }
        // cout << "full flow prefix: " << from_py_string(full_flow_prefix) << std::endl;
        wstring default_py_prefix(to_py_string(STR(FLOW123D_PYTHON_PREFIX)));
        // cout << "default py prefix: " << from_py_string(default_py_prefix) << std::endl;

        static wstring our_py_home(full_flow_prefix + ":" +default_py_prefix);
        Py_SetPythonHome( &(our_py_home[0]) );

#else
    (void)program_name; // not used
#endif //FLOW123D_PYTHON_PREFIX 

    // initialize the Python interpreter.
    py::initialize_interpreter();

    // update module path, first get current system path
    py::module_ sys = py::module_::import("sys");
    sys.attr("path").attr("append")(FLOW123D_SOURCE_DIR);
    std::string unit_tests_path = std::string(FLOW123D_SOURCE_DIR) + "/unit_tests";
    sys.attr("path").attr("append")( unit_tests_path.c_str() );
#ifdef FLOW123D_PYTHON_EXTRA_MODULES_PATH
    // than append flow123d Python modules path to sys.path
    std::string extra_paths(FLOW123D_PYTHON_EXTRA_MODULES_PATH);
    std::vector<std::string> extra_path_vec;
    boost::split(extra_path_vec, extra_paths, boost::is_any_of(":"));
    for (auto path : extra_path_vec)
        sys.attr("path").attr("append")(path.c_str());
#endif //FLOW123D_PYTHON_EXTRA_MODULES_PATH

    // call python and get paths available
//    //PyObject *moduleMain = PyImport_ImportModule("__main__");
//    py::module_ moduleMain = py::module_::import("__main__");
//    PyRun_SimpleString(python_sys_path.c_str());
//    //PyObject *func = PyObject_GetAttrString(moduleMain, "get_paths");
//    PyObject *func = moduleMain.attr("get_paths").cast<py::object>().release().ptr();
//    PyObject *args = PyTuple_New (0);
//    //PyObject *args = py::make_tuple(0, py::none(), "").cast<py::object>().release().ptr(); //py::tuple args
//    PyObject *result = PyObject_CallObject(func, args);
//    PythonLoader::check_error();
//
//    // save value so we dont have to call python again
//    PythonLoader::sys_path = string(PyUnicode_AsUTF8(result));
}



PythonRunning::~PythonRunning() {
    py::finalize_interpreter();
}

} // close namespace internal

#endif // FLOW123D_HAVE_PYTHON
