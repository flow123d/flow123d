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
        self.used_fields_dict = dict()
        self.result_fields_dict = dict()


    def __getattr__(self, attr):
        """ Allows direct access to items in 'used_fields_dict' field data dictionary.
            Example: use 'self.field_name' instead of 'self.used_fields_dict["field_name"]' """
        # print(attr)
        cache_data = self.used_fields_dict.get(attr, None)
        return cache_data[..., self.region_chunk_begin:self.region_chunk_end]
        

    def used_fields(self):
        """ Method allows to define list of fields used in evaluation. This is the default defintion
            that returns coords field and should be overwrite in descendant class. """
        field_list = ["X"]
        return field_list
        

    def repl(self, x):
        """ Method replicates scalar/vector/tensor field value to output vector. """
        return x[..., None]
        

    def _cache_reinit(self, data, result):
        """ Fill dictionary of input fields and result field """
        self.used_fields_dict.clear()
        for in_field in data:
            self.used_fields_dict[in_field.field_name()] = in_field.field_cache_array()
        
        self.result = result.field_name()
        self.result_fields_dict[self.result] = result.field_cache_array()


    def _cache_update(self, reg_chunk_begin, reg_chunk_end):
        """ Method called from cache_update in C++ code
            Needs to define __call__ method in descendant that executes evaluation """
        self.region_chunk_begin = reg_chunk_begin
        self.region_chunk_end = reg_chunk_end
        self.__call__()


    def _print_fields(self):
        """ Temporary method for development """
        print("Dictionary contains fields: ")
        print("1) Used fields: ")
        for key, val in self.used_fields_dict.items():
            print(key, ":")
            print(val)
        print("2) Result fields: ")
        for key, val in self.result_fields_dict.items():
            print(key, ":")
            print(val)

    