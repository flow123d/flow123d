#!/bin/python3
# author: David Flanderka

import sys
#import flowpy
import numpy as np

class PythonFieldBase():

    def __init__(self):
        """ Initialize object and its data members """
        self.region_chunk_begin = 0
        self.region_chunk_end = 300
        self.t = 0.0
        self.result = ""
        self.f_dict = dict()


    def __getattr__(self, attr):
        """ Allows direct access to items in 'f_dict' field data dictionary.
            Example: use 'self.field_name' instead of 'self.f_dict["field_name"]' """
        # print(attr)
        cache_data = self.f_dict.get(attr, None)
        return cache_data[:, self.region_chunk_begin:self.region_chunk_end]
        

    def _cache_reinit(self, data, result):
        """ Fill dictionary of input fields and result field """
        self.f_dict.clear()
        for in_field in data:
            self.f_dict[in_field.field_name()] = in_field.field_cache_array()
        
        self.result = result.field_name()
        self.f_dict[self.result] = result.field_cache_array()


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

    def _print_fields(self):
        """ Temporary method for development """
        print("Dictionary contains fields: ")
        for key, val in self.f_dict.items():
            print(key, ":")
            print(val)

    