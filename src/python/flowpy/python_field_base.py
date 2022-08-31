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

    # Method allows to define list of fields used in evaluation. This is the default defintion
    # that returns coords field and should be overwrite in descendant class. 
    def used_fields():
        field_list = ["X"]
        return field_list
        