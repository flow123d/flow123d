# This macro prints name and content of given variable.
#
# usage: PrintVar(FOO)


MACRO (PrintVar VAR)
    # put args into list
    message(STATUS "${VAR} =>" ${${VAR}})
ENDMACRO (PrintVar)