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

using namespace std;

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


PyObject * PythonLoader::load_module_from_file(const std::string& fname) {
    initialize();

    // don't know direct way how to load module form file, so we read it into string and use load_module_from_string
    std::ifstream file_stream;
    // check correct openning
    FilePath(fname, FilePath::input_file).open_stream(file_stream);
    file_stream.exceptions ( ifstream::failbit | ifstream::badbit );

    std::stringstream buffer;
    buffer << file_stream.rdbuf();

    string module_name;
    std::size_t pos = fname.rfind("/");
    if (pos != string::npos)
        module_name = fname.substr(pos+1);
    else
        module_name = fname;

    // cout << "python module: " << module_name <<endl;
    // DebugOut() << buffer.str() << "\n";
    // TODO: use exceptions and catch it here to produce shorter and more precise error message
    return load_module_from_string(module_name, buffer.str() );
}



PyObject * PythonLoader::load_module_from_string(const std::string& module_name, const std::string& source_string) {
    initialize();
    // compile string and check for errors
    PyObject * compiled_string = Py_CompileString(source_string.c_str(), module_name.c_str(), Py_file_input);
    PythonLoader::check_error();
    
    // import modle and check for error
    PyObject * result = PyImport_ExecCodeModule(module_name.c_str(), compiled_string);
    PythonLoader::check_error();
    return result;
}

PyObject * PythonLoader::load_module_by_name(const std::string& module_name) {
    initialize();

    // import module by dot separated path and its name
    PyObject * module_object = PyImport_ImportModule (module_name.c_str());
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



PyObject * PythonLoader::get_callable(PyObject *module, const std::string &func_name) {
    char func_char[func_name.size()+2];
    strcpy(func_char, func_name.c_str());
    PyObject * func = PyObject_GetAttrString(module, func_char );
    PythonLoader::check_error();

    if (! PyCallable_Check(func)) {
    	stringstream ss;
    	ss << "Field '" << func_name << "' from the python module: " << PyModule_GetName(module) << " is not callable." << endl;
    	THROW(ExcPythonError() << EI_PythonMessage( ss.str() ));
    }

    return func;
}



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
    Py_Initialize();

#ifdef FLOW123D_PYTHON_EXTRA_MODULES_PATH
    // update module path, first get current system path (Py_GetPath)
    // than append flow123d Python modules path to sys.path
    wstring path = Py_GetPath();
    path = path + to_py_string(":") + to_py_string(string(FLOW123D_PYTHON_EXTRA_MODULES_PATH));
    PySys_SetPath (path.c_str());
    
    
#endif //FLOW123D_PYTHON_EXTRA_MODULES_PATH

    // call python and get paths available
    PyObject *moduleMain = PyImport_ImportModule("__main__");
    PyRun_SimpleString(python_sys_path.c_str());
    PyObject *func = PyObject_GetAttrString(moduleMain, "get_paths");
    PyObject *args = PyTuple_New (0);
    PyObject *result = PyObject_CallObject(func, args);
    PythonLoader::check_error();
    
    // save value so we dont have to call python again
    PythonLoader::sys_path = string(PyUnicode_AsUTF8(result));
}



PythonRunning::~PythonRunning() {
    Py_Finalize();
}

} // close namespace internal

#endif // FLOW123D_HAVE_PYTHON
