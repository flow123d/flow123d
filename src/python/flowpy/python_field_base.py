#!/bin/python3
# author: David Flanderka

import sys

import flowpy

class PythonFieldBase(flowpy.PythonFieldBaseCPP):
    pass

    # Initialize object, set default values of region chunk range
    def __init__(self):
        flowpy.PythonFieldBaseCPP.__init__(self)
        self.region_chunk_begin = 0
        self.region_chunk_end = 300

    # Allows direct access to items in 'f_dict' field data dictionary.
    # Example: use 'self.field_name' instead of 'self.f_dict["field_name"]'
    def __getattr__(self, attr):
        cache_data = self.f_dict.get(attr, None)
        return cache_data[self.region_chunk_begin, self.region_chunk_end]
        
    # Allows call set_dict method from descendants
    def set_dict(self, data, result):
        super()._set_dict(data, result)

    # Method allows to define list of fields used in evaluation. This is the default defintion
    # that returns coords field and should be overwrite in descendant class. 
    def used_fields(self):
        field_list = ["X"]
        return field_list
        
    # Method called from cache_update in C++ code
    # Needs to define __call__ method in descendant that executes evaluation
    def cache_update(self, reg_chunk_begin, reg_chunk_end):
        self.region_chunk_begin = reg_chunk_begin
        self.region_chunk_end = reg_chunk_end
        self.__call__()
    