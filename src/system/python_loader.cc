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

#include "system/python_loader.hh"

//#include "system/system.hh"
#include "system/file_path.hh"
#include <string>
//#include <iostream>
#include <fstream>
//#include <sstream>
#include <pybind11.h>
#include <embed.h> // everything needed for embedding

using namespace std;
namespace py = pybind11;


void PythonLoader::initialize()
{
    static internal::PythonRunning _running;
}


py::module_ PythonLoader::load_module_from_file(const std::string& fname) {
    initialize();

    std::ifstream file_stream;
    FilePath f_path(fname, FilePath::input_file);
    f_path.open_stream(file_stream); // checks if file exists
    std::string parent_path = f_path.parent_path(); // add path to PythonPath
    PythonLoader::add_sys_path(parent_path);

    string module_name = f_path.stem();
    py::module_ module;
    try {
    	module = py::module_::import(module_name.c_str());
    } catch (const py::error_already_set &ex) {
        PythonLoader::throw_error(ex);
    }
    return module;
}



py::module_ PythonLoader::load_module_from_string(const std::string& module_name, const std::string& func_name, const std::string& source_string) {
    initialize();

    py::module_ user_module;
    py::dict globals = py::globals();

    try {
        py::exec(source_string.c_str(), globals, globals);
        user_module = py::module::create_extension_module(module_name.c_str(), "", new PyModuleDef());
        user_module.add_object(func_name.c_str(), globals[func_name.c_str()]);
    } catch (const py::error_already_set &ex) {
        PythonLoader::throw_error(ex);
    }
    return user_module;
}

py::module_ PythonLoader::load_module_by_name(const std::string& module_name) {
    initialize();

    // import module by dot separated path and its name
    py::module_ module_object;
    try {
        module_object = py::module_::import (module_name.c_str());
    } catch (const py::error_already_set &ex) {
        PythonLoader::throw_error(ex);
    }

    return module_object;
}


void PythonLoader::throw_error(const py::error_already_set &ex) {
    PyObject *ptype = ex.type().ptr();
    PyObject *pvalue = ex.value().ptr();
    PyObject *ptraceback = ex.trace().ptr();

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
    auto path_vec = PythonLoader::get_python_path();
    string python_path;
    for (auto p : path_vec) {
        python_path += p + "\n";
    }

    // construct error message
    string py_message =
               "\nType: " + string(PyUnicode_AsUTF8(PyObject_Str( ex.type().ptr() ))) + "\n"
             + "Message: " + string(PyUnicode_AsUTF8(PyObject_Str( ex.value().attr("args").ptr() ))) + "\n"
             + "Traceback: \n" + str_traceback + "\n"
             + "Paths: " + "\n" + python_path + "\n";
    
    THROW(ExcPythonError() << EI_PythonMessage( py_message ));
}



void PythonLoader::add_sys_path(const std::string &path)
{
    py::module_ sys = py::module_::import("sys");
    auto sys_paths = sys.attr("path");
    // Checks if path exists
    for (auto sp : sys_paths) {
        std::string sys_path = sp.cast<std::string>();
        if (sys_path == path) {
            WarningOut() << "Path '" << sys_path << "' already exists in Python sys path and cannot be added repeatedly!";
            return;
        }
    }

    // Adds only one time
    sys_paths.attr("append")(path);
}

std::vector<std::string> PythonLoader::get_python_path()
{
    py::module_ sys = py::module_::import("sys");
    auto sys_paths = sys.attr("path");
    std::vector<std::string> path_vec;
    for (auto sp : sys_paths) {
        std::string sys_path = sp.cast<std::string>();
        path_vec.push_back(sys_path);
    }

    return path_vec;
}



namespace internal {

PythonRunning::PythonRunning()
{
    // initialize the Python interpreter.
    py::initialize_interpreter();
    std::string flowpy_path = std::string(FLOW123D_SOURCE_DIR) + "/src/python/flowpy";
    std::string fieldproxypy_path = std::string(FLOW123D_SOURCE_DIR) + "/build_tree/src";
    PythonLoader::add_sys_path(flowpy_path);
    PythonLoader::add_sys_path(fieldproxypy_path);

#ifdef FLOW123D_PYTHON_EXTRA_MODULES_PATH
    // update module path, append flow123d Python modules path to sys.path
    std::stringstream extra_paths(FLOW123D_PYTHON_EXTRA_MODULES_PATH);
    std::string extra_path;
    while ( std::getline(extra_paths, extra_path, ':') )
    {
       PythonLoader::add_sys_path(extra_path);
    }
#endif //FLOW123D_PYTHON_EXTRA_MODULES_PATH
}




PythonRunning::~PythonRunning() {
    py::finalize_interpreter();
}

} // close namespace internal
