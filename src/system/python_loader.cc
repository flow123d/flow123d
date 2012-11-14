/*
 * python_loader.cc
 *
 *  Created on: Aug 31, 2012
 *      Author: jb
 */


#ifdef HAVE_PYTHON

#include "system/python_loader.hh"
#include "global_defs.h"
#include "system/system.hh"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

PyObject * PythonLoader::load_module_from_file(const std::string& fname) {
    // don't know direct way how to load module form file, so we read it into string and use load_module_from_string
    std::ifstream file_stream( fname.c_str() );
    std::stringstream buffer;
    buffer << file_stream.rdbuf();

    string module_name;
    unsigned int pos = fname.rfind("/");
    if (pos != string::npos)
        module_name = fname.substr(pos+1);
    else
        module_name = fname;

    cout << "python module: " << module_name <<endl;
    // TODO: use exceptions and catch it here to produce shorter and more precise error message
    return load_module_from_string(module_name, buffer.str() );
}



PyObject * PythonLoader::load_module_from_string(const std::string& module_name, const std::string& source_string) {

    // for unknown reason module name is non-const, so we have to make char * copy
    char * tmp_name = new char[ module_name.size() + 2 ];
    strcpy( tmp_name, module_name.c_str() );
    PyObject * result = PyImport_ExecCodeModule(tmp_name,
            Py_CompileString( source_string.c_str(), "flow123d_python_loader", Py_file_input ) );

    if (result == NULL) {
        PyErr_Print();
        std::cerr << "Error: Can not load python module '" << module_name << "' from string:" << std::endl;
        std::cerr << source_string << std::endl;
    }

    delete[] tmp_name;
    return result;
}


internal::PythonRunning PythonLoader::py_running;


namespace internal {

PythonRunning::PythonRunning() {
    DBGMSG("Init python\n");
    Py_Initialize();
    DBGMSG("Init python, done\n");
}



PythonRunning::~PythonRunning() {
    Py_Finalize();
}

} // close namespace internal

#endif // HAVE_PYTHON
