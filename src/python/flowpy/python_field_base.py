#!/bin/python3
# author: David Flanderka

import types
from typing import *
import numpy as np
from flow123d_python_api import FieldCacheProxy     # neccessary for conversion to numpy array

class PythonFieldBase:
    """
    Base class for user field evaluation classes. The user class defines the evaluation method
    for every field to evaluate. This evaluation can depend on other fields specified by the key `used_fields` of
    the `FieldPython` input record of Flow123d. The evaluation through cal of `_cache_update` evaluates many quadrature
    points (n_points) at once. The field values are stored in the numpy arrays where the last dimension
    corresponds to the quadrature points.
    E.g.  scalar field values are stored in 1D array of shape (n_points)
    3D vector field values are stored in array of shape (3, n_points)
    tensor values are stored in array of shape (3, 3, n_points)



    Usage:
    YAMLFile:
        water_source_density: !FieldPython
          source_file: eval.py
          class: Eval
          used_fields: ["X", "cross_section"]
        conductivity: !FieldPython
          source_file: eval.py
          class: Eval
          used_fields: ["X"]

    eval.py:

    class UserEvaluation(PythonFieldBase):
        # User class will be instantiated just once and shared for among all evaluated fields.
        # That allows common precalculations.
        # PythonFieldBase also allows dot acces to the input fields.
        def water_source_density(self):
            # evaluation of the `field_b` using self as dictionary of the input fields
            # specified in the YAML file
            return np.norm(self.X) * self.cross_section
        
        def conductivity(self):
            # evaluation of
            z = self.X[2]
            return np.exp(-z)
    """
    
    # Singleton like instances of the user field evaluation classes.
    _instances = dict()

    @staticmethod
    def _create(module: types.ModuleType, class_name: str) -> 'PythonFieldBase':
        """ 
        Create and return unique instance of the user field evaluation class of name `class_name` from `module`.
        """
        full_name = f"{module.__name__}.{class_name}"
        try:
            instance = PythonFieldBase._instances[full_name]
        except KeyError:
            class_ = getattr(module, class_name)
            instance = class_()
            PythonFieldBase._instances[full_name] = instance
        return instance

    def __init__(self):
        """ 
        Constructor. Define all generic attributes.
        """
        # Slice of the quadrature points to evaluate, set by _cache_update.
        self._region_chunk_begin = 0
        self._region_chunk_end = 0
        # Currently evaluated time.
        self.t = 0.0
        # Dictionary of the input fields. Access through dot syntax:
        # self.<input_field>
        self._used_fields_dict = dict()
        # Dictionary of the output fields. No direct access. Reuslts of the user methods are stored
        # in `_cache_update`.
        self._result_fields_dict = dict()

    def __getattr__(self, attr):
        """
        Access to the input fields through th e"dot" sysntax.
        Example: use 'self.field_name' id equal to 'self.used_fields_dict["field_name"]'
        """
        cache_data = self._used_fields_dict.get(attr, None)
        return cache_data[..., self._region_chunk_begin:self._region_chunk_end]
        
    @staticmethod
    def repl(x):
        """
        Replicates the input array along the axis of the quadrature points (the last one).
        Is necessary for the scalar arrays, e.g.
        ```
             y = self.X[1]
             return self.repl( np.eye(3) ) * y
        ```
        Without repl, the expression `np.eye(3) * y` would fail as we try to multiply constant tensor of shape (3,3)
        by the scalar field of the shape (n_points). Using `repl` we multiply shape (3,3,1) by shape
        (n_points) which is broadcasted to common the shape (3,3,n_points).
        """
        return x[..., None]
        

    def _cache_reinit(self, time: float, data: List[FieldCacheProxy], result: FieldCacheProxy) -> None:
        """
        Create arrays as wrappers to given C++ field value caches passed as FieldCacheProxy.
        One reinit is called for every result field.
        """
        self._used_fields_dict.clear()
        for in_field in data:
            in_array = np.array(in_field, copy=False)
            in_array.flags.writeable = False
            self._used_fields_dict[in_field.field_name()] = in_array
        self._result_fields_dict[result.field_name()] = np.array(result, copy=False)
        
        self.t = time


    def _cache_update(self, field_name: str, reg_chunk_begin: int, reg_chunk_end: int):
        """
        Method called from cache_update in C++ code.
        Needs to define the method with same name as name of evaluated field in descendant
        that executes evaluation.
        """
        self._region_chunk_begin = reg_chunk_begin
        self._region_chunk_end = reg_chunk_end
        res_array = getattr(self, field_name)()
        
        # Check number of dimensions and shape of result array
        result_shape = res_array.shape
        expect_shape = self._result_fields_dict[field_name][..., self._region_chunk_begin:self._region_chunk_end].shape
        if result_shape != expect_shape:
            raise ValueError(f"Invalid shape of '{field_name}' method result. Must be {expect_shape}.")

        self._result_fields_dict[field_name][..., self._region_chunk_begin:self._region_chunk_end] = res_array


    def _print_fields(self):
        """ Auxiliary method for development """
        print("Dictionary contains fields: ")
        print("1) Used fields: ")
        for key, val in self._used_fields_dict.items():
            print(key, ":")
            print(val)
        print("2) Result fields: ")
        for key, val in self._result_fields_dict.items():
            print(key, ":")
            print(val)

    
