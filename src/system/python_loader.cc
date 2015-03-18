/*
 * python_loader.cc
 *
 *  Created on: Aug 31, 2012
 *      Author: jb
 */

#include "global_defs.h"

#ifdef FLOW123D_HAVE_PYTHON

#include "system/python_loader.hh"

#include "system/system.hh"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

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

    PyObject * compiled_string = Py_CompileString( source_string.c_str(), module_name.c_str(), Py_file_input );
    PythonLoader::check_error();

    PyObject * result = PyImport_ExecCodeModule(tmp_name, compiled_string);
    PythonLoader::check_error();

    delete[] tmp_name;
    return result;
}



void PythonLoader::check_error() {
    if (PyErr_Occurred()) {
    	PyObject *ptype, *pvalue, *ptraceback, *pystr;
		PyErr_Fetch(&ptype, &pvalue, &ptraceback);

		if ( PyString_Check(pvalue) ) { // error message contains simple string
			THROW(ExcPythonError() << EI_PythonMessage( PyString_AsString(pvalue) ));
		} else { // set of string values
			stringstream ss;
			PyObject* err_type = PySequence_GetItem(pvalue, 0);
			pystr = PyObject_Str(err_type);
			ss << PyString_AsString(pystr) << endl;
			PyObject* descr = PySequence_GetItem(pvalue, 1);
			PyObject* file = PySequence_GetItem(descr, 0);
			pystr = PyObject_Str(file);
			ss << "  File \"" << PyString_AsString(pystr) << "\"";
			PyObject* row = PySequence_GetItem(descr, 1);
			pystr = PyObject_Str(row);
			ss << ", line " << PyString_AsString(pystr) << endl;
			PyObject* line = PySequence_GetItem(descr, 3);
			pystr = PyObject_Str(line);
			ss << PyString_AsString(pystr);
			PyObject* col = PySequence_GetItem(descr, 2);
			pystr = PyObject_Str(col);
			int pos = atoi( PyString_AsString(pystr) );
			ss << setw(pos) << "^" << endl;
			THROW(ExcPythonError() << EI_PythonMessage( ss.str() ));
		}
    }
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

// currently we support only Python 2.7
//
#if FLOW123D_PYTHONLIBS_VERSION_MAJOR<3
    #define to_py_string      string
    #define from_py_string    string
    #define PY_STRING string
#else
    #define PY_STRING wstring
#endif

#define STR_EXPAND(tok) #tok
#define STR(tok) string(STR_EXPAND(tok))

namespace internal {

PythonRunning::PythonRunning(const std::string& program_name)
{
#ifdef FLOW123D_PYTHON_PREFIX
        static PY_STRING _python_program_name = to_py_string(program_name);
        Py_SetProgramName( &(_python_program_name[0]) );
        PY_STRING full_program_name = Py_GetProgramFullPath();
        cout << "full program name: " << from_py_string(full_program_name) << std::endl;

        size_t pos = full_program_name.rfind( to_py_string("flow123d") );
        DBGMSG("pos: %d\n", pos);
        ASSERT(pos != PY_STRING::npos, "non flow123d binary");
        PY_STRING full_flow_prefix=full_program_name.substr(0,pos-string("/bin/").size() );
        cout << "full flow prefix: " << from_py_string(full_flow_prefix) << std::endl;
        PY_STRING default_py_prefix(to_py_string(STR(FLOW123D_PYTHON_PREFIX)));
        cout << "default py prefix: " << from_py_string(default_py_prefix) << std::endl;

        static PY_STRING our_py_home(full_flow_prefix + ":" +default_py_prefix);
        Py_SetPythonHome( &(our_py_home[0]) );
        
        /*
        Py_GetPath();

        static PY_STRING our_py_path;
        string python_subdir("/lib/python" + STR(FLOW123D_PYTHONLIBS_VERSION_MAJOR) + "." + STR(FLOW123D_PYTHONLIBS_VERSION_MINOR));
        our_py_path+=full_flow_prefix + to_py_string( python_subdir + "/:");
        our_py_path+=full_flow_prefix + to_py_string( python_subdir + "/plat-cygwin:");
        our_py_path+=full_flow_prefix + to_py_string( python_subdir + "/lib-dynload:");
        our_py_path+=default_py_prefix + to_py_string( python_subdir + "/lib-dynload:");
        our_py_path+=default_py_prefix + to_py_string( python_subdir + "/lib-dynload:");
        our_py_path+=default_py_prefix + to_py_string( python_subdir + "/lib-dynload:");

        Py_SetPath(our_py_path);
        //our_py_path+=full_flow_prefix + to_py_string( python_subdir);

//        string prefix = ;
        */
        cout << "Python path: " << from_py_string( Py_GetPath() ) << std::endl;
        cout << "Python home: " << from_py_string( Py_GetPythonHome() ) << std::endl;
        cout << "Python prefix: " << from_py_string( Py_GetPrefix() ) << std::endl;
        cout << "Python exec prefix: " << from_py_string( Py_GetExecPrefix() ) << std::endl;
        
        // 1. set program name
        // 2. get prefix
        // 3. get python full path
        // 4. extract prefix from the full path
        // 5. substitute the prefix in python path
        // 6. set python path (suppress its computation during Py_Initialize)
        /*
        wchar_t wbuf[ python_home.size() ];
        size_t num_chars = mbstowcs( wbuf, python_home.c_str(), python_home.size() );
	static wstring _python_home_storage( wbuf, num_chars ); // permanent non-const storage required

        std::wcout << "new python home: " << _python_home_storage << std::endl;

	Py_SetProgramName( &(_python_home_storage[0]) );

        char buff[ 1024 ];
        num_chars = wcstombs(buff, Py_GetPath(), 256);
        buff[1024]=0;
        string str(buff, num_chars);
        std::cout << "Python path: " << str << std::endl;



                wchar_t wbuf[ python_home.size() ];
        size_t num_chars = mbstowcs( wbuf, python_home.c_str(), python_home.size() );
        static wstring _python_home_storage( wbuf, num_chars ); // permanent non-const storage required

        std::wcout << "new python home: " << _python_home_storage << std::endl;

        Py_SetProgramName( &(_python_home_storage[0]) );

        char buff[ 1024 ];
        num_chars = wcstombs(buff, Py_GetPath(), 1024);
        //std::cout << "num chars: " << num_chars << std::endl;
        std::cout << "Python path: " << buff << std::endl;
        //num_chars = wcstombs(buff, Py_GetPythonHome(), 1024);
        //std::cout << "Python home: " << buff << std::endl;
        num_chars = wcstombs(buff, Py_GetProgramName(), 1024);
        std::cout << "Python prog. name: " << buff << std::endl;
        num_chars = wcstombs(buff, Py_GetPrefix(), 1024);
        std::cout << "Python prefix: " << buff << std::endl;
        num_chars = wcstombs(buff, Py_GetProgramFullPath(), 1024);
        std::cout << "Python full: " << buff << std::endl;
*/
#endif
    Py_Initialize();
}



PythonRunning::~PythonRunning() {
    Py_Finalize();
}

} // close namespace internal

#endif // FLOW123D_HAVE_PYTHON
