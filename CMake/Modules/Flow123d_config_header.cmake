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


# macro which replaces add_definitions function and stores ONLY numeric constants in the DEFINITIONS_LIST
# also generates DEFINITIONS_CONTENT variable containing valid header file definitions
# MACRO also adds prefix FLOW123D_
# usage: 
# FLOW_DEFINE (PYTHONLIBS_VERSION_MAJOR 3)      -> #define FLOW123D_PYTHONLIBS_VERSION_MAJOR 3
MACRO (FLOW_DEFINE_CONSTANT variable_name numeric_value)   
    MESSAGE (STATUS "Adding: '${variable_name}' ${numeric_value}")
    LIST (APPEND DEFINITIONS_LIST "FLOW123D_${variable_name}=${numeric_value}")
    SET (DEFINITIONS_CONTENT "${DEFINITIONS_CONTENT}#define FLOW123D_${variable_name} ${numeric_value}\n")
ENDMACRO (FLOW_DEFINE_CONSTANT variable_name numeric_value)



# MACRO will generate definitions.tmp file in which all definition will be stored
# after that python script will expand this file into valid header file config.h
MACRO (GENERATE_CONFIG_H file_path)
    # prepare tmp output location variable 
    MESSAGE ("GENERATING CONFIG_H_TMP")
    set(OUTPUT_TMP_PATH "${file_path}.tmp")
    # write definitions to tmp file
    FILE (WRITE "${OUTPUT_TMP_PATH}" "${DEFINITIONS_CONTENT}")
    execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different ${OUTPUT_TMP_PATH} ${file_path})
    execute_process(COMMAND ${CMAKE_COMMAND} -E remove ${OUTPUT_TMP_PATH})
ENDMACRO (GENERATE_CONFIG_H file_path)