#!/bin/python3
# author: David Flanderka

import sys
import fieldproxypy
import numpy as np
from typing import *


class PythonFieldBase():
    _instances = dict()

    @staticmethod
    def _create(module, class_name):
        """ Creates instance of class_name if doesn't exist. Stores its to _instances and returns. """
        if class_name not in PythonFieldBase._instances:
            # module = __import__(module_name)
            class_ = getattr(module, class_name)
            PythonFieldBase._instances[class_name] = class_()
        return PythonFieldBase._instances.get(class_name, None)


    def __init__(self):
        """ Initialize object and its data members """
        self.region_chunk_begin = 0
        self.region_chunk_end = 300
        self.t = 0.0
        self.used_fields_dict = dict()
        self.result_fields_dict = dict()


    def __getattr__(self, attr):
        """ Allows direct access to items in 'used_fields_dict' field data dictionary.
            Example: use 'self.field_name' instead of 'self.used_fields_dict["field_name"]' """
        cache_data = self.used_fields_dict.get(attr, None)
        return cache_data[..., self.region_chunk_begin:self.region_chunk_end]
        

    def repl(self, x):
        """ Method replicates x value to size of evaluated vector. 
            Example, use this method for fill result:
             > y = self.X[1]
             > return self.repl( np.eye(3) ) * y
            In this example, method creates vector of identity matrices. Size of vector 
            is same as size of field result. This fact is ensured by this method. Then 
            vector of identity matrices id multiplicated by values on y-coord.
            """
        return x[..., None]
        

    def _cache_reinit(self, time: float, data: List['FieldCacheProxy'], result: 'FieldCacheProxy') -> None:
        """ Fill dictionary of input fields and result field, set time """
        self.used_fields_dict.clear()
        for in_field in data:
            in_array = np.array(in_field, copy=False)
            in_array.flags.writeable = False
            self.used_fields_dict[in_field.field_name()] = in_array
        self.result_fields_dict[result.field_name()] = np.array(result, copy=False)
        
        self.t = time


    def _cache_update(self, field_name: str, reg_chunk_begin: int, reg_chunk_end: int):
        """ Method called from cache_update in C++ code.
            Needs to define the method with same name as name of evaluated field in descendant 
            that executes evaluation. """
        self.region_chunk_begin = reg_chunk_begin
        self.region_chunk_end = reg_chunk_end
        res_array = getattr(self, field_name)()
        
        # Check number of dimensions and shape of result array
        result_shape = res_array.shape
        expect_shape = self.result_fields_dict[field_name][..., self.region_chunk_begin:self.region_chunk_end].shape
        if result_shape!=expect_shape:
            raise ValueError(f"Invalid shape of '{field_name}' method result. Must be {expect_shape}.")
        
        self.result_fields_dict[field_name][..., self.region_chunk_begin:self.region_chunk_end] = res_array


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

    