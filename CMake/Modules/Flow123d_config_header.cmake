# variable for storing definitions
SET(DEFINITIONS_LIST "")
SET(DEFINITIONS_CONTENT "")

# macro which replaces add_definitions function and store all values to the DEFINITIONS_LIST
# also generates DEFINITIONS_CONTENT variable containing valid header file definitions
# MACRO also adds prefix FLOW123D_
# usage: 
# FLOW_DEFINE (Flow123d_DEBUG)      -> #define FLOW123D_Flow123d_DEBUG 1
# FLOW_DEFINE (Foo)                 -> #define FLOW123D_Foo 1
# FLOW_DEFINE (Bar "string value")  -> #define FLOW123D_Bar "string value"
MACRO (FLOW_DEFINE variable_name)	
    SET(variable_value "${ARGN}")

    if ("${variable_value}" MATCHES "^$")
    	MESSAGE (STATUS "Adding: '${variable_name}' 1")
    	LIST (APPEND DEFINITIONS_LIST "FLOW123D_${variable_name}=1")
        SET (DEFINITIONS_CONTENT "${DEFINITIONS_CONTENT}#define FLOW123D_${variable_name} 1\n")
    else()
    	MESSAGE (STATUS "Adding: '${variable_name}' '${variable_value}'")
		LIST (APPEND DEFINITIONS_LIST "FLOW123D_${variable_name}=${variable_value}")
        SET (DEFINITIONS_CONTENT "${DEFINITIONS_CONTENT}#define FLOW123D_${variable_name} \"${variable_value}\"\n")
    endif()
ENDMACRO (FLOW_DEFINE variable_name)


# MACRO will generate definitions.tmp file in which all definition will be stored
# after that python script will expand this file into valid header file config.h
MACRO (GENERATE_CONFIG_H file_path)
    # write definitions to tmp file
    MESSAGE ("GENERATING CONFIG_H_TMP")
    FILE (WRITE "${file_path}" "${DEFINITIONS_CONTENT}")
ENDMACRO (GENERATE_CONFIG_H)

# SENSITIVE DEFINITIONS
# FLOW123D_HAVE_PETSC:
#       CMakeLists.txt
#       application_base.cc         src/system
#       application_base.hh         src/system
# FLOW123D_HAVE_MPI:
#       CMakeLists.txt
# ARMA_NO_DEBUG
#       
# FLOW123D_HAVE_CXX11_FULL:
#       CMakeLists.txt
#       system.hh                   src/system
#       sys_profiler.hh             src/system
# FLOW123D_HAVE_CXX11_DRAFT:
#       CMakeLists.txt
#       sys_profiler.hh             src/system


# Other definitions
# FLOW123D_DEBUG_PROFILER 
#       CMakeLists.txt
#       config.cmake.template
#       global_defs.h               src/system
#       profiler_test.cpp           unit_tests/system
#       sys_profiler.cc             src/system
#       sys_profiler.hh             src/system
# FLOW123D_RUN_UNIT_BENCHMARKS
#       CMakeLists.txt
#       config.cmake.template 
#       field_speed_test.cpp        unit_tests/fields
#       region_test.cpp             unit_tests/mesh
#       tokenizer_speed_test.cpp    unit_tests/system
# FLOW123D_HAVE_SINCOS
#       CMakeLists.txt
#       fpaux.hh                    third_party/fparser-4.5.1/extrasrc
# FLOW123D_CYGWIN
#       CMakeLists.txt
#       equation.cc                 src/coupling
#       file_path.cc                src/system
# HAVE_EXEC_INFO
#       CMakeLists.txt
#       exceptions.cc               src/system
# HAVE_DEMAGLER
#       exceptions.cc               src/system   
# HAVE_PYTHON
#       CMakeLists.txt
#       configure                   third_party/gtest-1.7.0
#       configure.ac                third_party/gtest-1.7.0
#       field_algo_base.impl.hh     src/fields
#       field_python.hh             src/fields
#       field_python.impl.hh        src/fields
#       field_python_test.cpp       unit_tests/fields
#       field_speed_test.cpp        unit_tests/fields
#       main.cc                     src
#       Makefile.am                 third_party/gtest-1.7.0
#       python_loader.cc            src/system
#       python_loader.hh            src/system
# PYTHONLIBS_VERSION_...
#       CMakeLists.txt
#       python_loader.cc            src/system
# PYTHON_PREFIX
#       CMakeLists.txt
#       python_loader.cc            src/system


