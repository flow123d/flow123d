# This macro sets first applicable/valid value passed in call
#
# SET_VALID_VALUE (VARIABLE_NAME [OPTION_1 [OPTION_2 [ OPTION_3 ...]]]
# variable with name VARIABLE_NAME will have value of OPTION_1,
# if OPTION_1 value is empty then OPTION_2 and so on

MACRO (SET_VALID_VALUE)
    # put args into list
    SET(tmp_list_var "${ARGN}")
    LIST (GET tmp_list_var 0 tmp_var_name)
    LIST (REMOVE_AT tmp_list_var 0)

    # single-step iteration to get first non empty value
    FOREACH(loop_var IN LISTS tmp_list_var)
        SET("${tmp_var_name}" "${loop_var}")
        BREAK()
    ENDFOREACH()
ENDMACRO (SET_VALID_VALUE)