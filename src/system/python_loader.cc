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

void PythonLoader::initialize(const std::string &python_home)
{
    static internal::PythonRunning _running(python_home);
}


PyObject * PythonLoader::load_module_from_file(const std::string& fname) {
    initialize();

    // don't know direct way how to load module form file, so we read it into string and use load_module_from_string
    std::ifstream file_stream( fname.c_str() );
    // check correct openning
    INPUT_CHECK(! file_stream.fail(), "Can not open input file '%s'.\n", fname.c_str() );
    file_stream.exceptions ( ifstream::failbit | ifstream::badbit );

    std::stringstream buffer;
    buffer << file_stream.rdbuf();

    string module_name;
    unsigned int pos = fname.rfind("/");
    if (pos != string::npos)
        module_name = fname.substr(pos+1);
    else
        module_name = fname;

    // cout << "python module: " << module_name <<endl;
    // DBGMSG("%s\n", buffer.str().c_str());
    // TODO: use exceptions and catch it here to produce shorter and more precise error message
    return load_module_from_string(module_name, buffer.str() );
}



PyObject * PythonLoader::load_module_from_string(const std::string& module_name, const std::string& source_string) {
    initialize();

    // for unknown reason module name is non-const, so we have to make char * copy
    char * tmp_name = new char[ module_name.size() + 2 ];
    strcpy( tmp_name, module_name.c_str() );
    PyObject * compiled_string = Py_CompileString( source_string.c_str(), "flow123d_python_loader", Py_file_input );
    if (! compiled_string) {
        PyErr_Print();
        std::cerr << "Error: Can not compile python string:\n'" << source_string << std::endl;
    }
    PyObject * result = PyImport_ExecCodeModule(tmp_name, compiled_string);

    if (result == NULL) {
        PyErr_Print();
        std::cerr << "Error: Can not load python module '" << module_name << "' from string:" << std::endl;
        std::cerr << source_string << std::endl;
    }

    delete[] tmp_name;
    return result;
}



namespace internal {

PythonRunning::PythonRunning(const std::string& python_home)
{
    /*
     * Possibly set alternative path to non-library part of Python.
     * This is necessary for release builds.
     */
    if (python_home != "") {
#if PYTHON_VERSION<3
		static string _python_home_storage = python_home; // permanent non-const storage required
		//Py_SetPythonHome( &(_python_home_storage[0]) );
		//Py_SetPythonHome( ".." );
		Py_SetProgramName( &(_python_home_storage[0]) );
#else
    	wchar_t *buf = new wchar_t[ python_home.size() ];
		size_t num_chars = mbstowcs( buf, python_home.c_str(), python_home.size() );
		wstring ws( buf, num_chars );
		delete[] buf;
		Py_SetProgramName( &(ws[0]) );

#endif
        std::cout << Py_GetPath() << std::endl;
    }
    Py_Initialize();
}



PythonRunning::~PythonRunning() {
    Py_Finalize();
}

} // close namespace internal

#endif // HAVE_PYTHON
