#!/bin/python3
# author: David Flanderka

import sys

import flowpy

class PythonFieldBase(flowpy.PythonFieldBaseCPP):
    pass

    # Allows direct access to items in 'f_dict' field data dictionary.
    # Example: use 'self.field_name' instead of 'self.f_dict["field_name"]'
    def __getattr__(self, attr):
        return self.f_dict.get(attr, None)
        