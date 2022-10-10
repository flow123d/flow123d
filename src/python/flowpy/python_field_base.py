#!/bin/python3
# author: David Flanderka

import sys

import flowpy

class PythonFieldBase(flowpy.PythonFieldBaseCPP):
    pass

    def __init__(self):
        """ Initialize object, set default values of region chunk range """
        flowpy.PythonFieldBaseCPP.__init__(self)
        self.region_chunk_begin = 0
        self.region_chunk_end = 300

    def __getattr__(self, attr):
        """ Allows direct access to items in 'f_dict' field data dictionary.
            Example: use 'self.field_name' instead of 'self.f_dict["field_name"]' """
        cache_data = self.f_dict.get(attr, None)
        return cache_data[:, self.region_chunk_begin:self.region_chunk_end]
        
    def set_dict(self, data, result):
        """ Allows call set_dict method from descendants """
        super()._set_dict(data, result)

    def used_fields(self):
        """ Method allows to define list of fields used in evaluation. This is the default defintion
            that returns coords field and should be overwrite in descendant class. """
        field_list = ["X"]
        return field_list
        
    def cache_update(self, reg_chunk_begin, reg_chunk_end):
        """ Method called from cache_update in C++ code
            Needs to define __call__ method in descendant that executes evaluation """
        self.region_chunk_begin = reg_chunk_begin
        self.region_chunk_end = reg_chunk_end
        self.__call__()
    