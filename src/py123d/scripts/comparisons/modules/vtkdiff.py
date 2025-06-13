#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Brezina
# ----------------------------------------------

import sys
import os 
import math
import vtk
import vtk.util.numpy_support as nps
import numpy as np
import argparse

from py123d.scripts.comparisons.modules import InPlaceComparison
from py123d.scripts.core.base import Paths, IO


def get_parsed_args():
    parser = argparse.ArgumentParser(
        prog="vtk_diff",
        description="Numerical comparison of two VTK/VTU files: reference vs. test result.")
    parser.add_argument("-a", "--abs_tol", dest='atol', type=float, default=1e-5, help="Absolute tolerance.")
    parser.add_argument("-r", "--rel_tol", dest='rtol', type=float, default=1e-2, help="Relative tolerance.")
    parser.add_argument("-i", "--interpolate", type=bool, default=False, help="Allow different meshes, interpolate to the reference mesh first.")
    parser.add_argument('ref_file',
                        help="Reference VTK/VTU file.")
    parser.add_argument('test_file',
                        help="Test VTK/VTU file.")

    a = parser.parse_args()
    print(a)
    return a

class VTKDiff(InPlaceComparison):

    def error(self, message):
        self.output.write(message)
        self.have_match = False
        
    def check_nans(dataset):
        names = self.get_array_names(dataset)
        for name in names:
            array = nps.vtk_to_numpy(ref_data.GetArray(names.index(name)))
            if np.any(np.isnan(array)):
                self.error(f"Reference array {name} contains NaNs.") 
                return True
                
        return False

    def compare_pipeline(self, f_ref, f_test):
        vtu_ref = self.read_vtu(f_ref)
        if check_nans(vtu_ref):
            return 
        vtu_test = self.read_vtu(f_test)
        if check_nans(vtu_test):
            return 
        
        if self.interpolate:
            vtu_ref.Update()
            vtu_test.Update()
            vtu_ref, vtu_test = self.interpolate_filter(vtu_ref, vtu_test)

        vtu_ref.Update()
        vtu_test.Update()
        #print(vtu_ref.GetOutput())
        #print(vtu_test.GetOutput())
        self.output.write('Cell data:')
        self.compare_data(vtu_ref.GetOutput().GetCellData(), vtu_test.GetOutput().GetCellData())
        self.output.write('Point data:')
        self.compare_data(vtu_ref.GetOutput().GetPointData(), vtu_test.GetOutput().GetPointData())

    def read_vtu(self, fname):
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(fname)
        return reader

    def clean_grid(self, vtk_data):
        point_data = vtk.vtkCellDataToPointData()
        point_data.SetInputData(vtk_data.GetOutput())
        point_data.Update()

        # Unfortunately this filter is able to compute average of merged points
        # This cause spurious differences for permuted meshes.
        clean_grid = vtk.vtkExtractUnstructuredGrid()
        clean_grid.SetInputData(point_data.GetOutput())
        clean_grid.ExtentClippingOff()
        clean_grid.CellClippingOff()
        clean_grid.PointClippingOff()
        clean_grid.MergingOn()
        clean_grid.Update()

        return clean_grid


    def interpolate_filter(self, ref, test):
        ref_clean = self.clean_grid(ref)
        test_clean = self.clean_grid(test)

        #print(ref.GetOutput().GetNumberOfPoints())
        #print(ref_point_data.GetOutput().GetNumberOfPoints())
        #print(ref_point_data.GetValidPoints().GetNumberOfTuples());

        test_on_ref = vtk.vtkProbeFilter()
        test_on_ref.SetInputData(ref_clean.GetOutput()) # points
        test_on_ref.SetSourceData(test_clean.GetOutput()) # data
        test_on_ref.Update();
        #print(test_on_ref.GetOutput().GetNumberOfPoints())
        #print(test_on_ref.GetValidPoints().GetNumberOfTuples());
        return ref_clean, test_on_ref

    def compare_data(self, ref_data, test_data):
        ref_names = self.get_array_names(ref_data)
        #print(kind, "ref arrays:", ref_names)
        test_names = self.get_array_names(test_data)
        try:
            test_names.remove('vtkValidPointMask')
        except Exception:
            pass
        #print(kind, "test arrays:", test_names)
        isec_names = set(ref_names).intersection(set(test_names))
        if len(isec_names) < len(ref_names):
            self.error(f"DIFF Missing fields in test: {set(ref_names).difference(set(test_names))}")
        if len(isec_names) < len(test_names):
            self.error(f"DIFF Surplus fields in test: {set(test_names).difference(set(isec_names))}")
        for name in isec_names:
            #print(name)
            ref_array = nps.vtk_to_numpy(ref_data.GetArray(ref_names.index(name)))
            test_array = nps.vtk_to_numpy(test_data.GetArray(test_names.index(name)))
            self.diff_arrays(ref_array, test_array, name)

    def get_array_names(self, data):
        return [data.GetArrayName(i) for i in range(data.GetNumberOfArrays())]

    def diff_arrays(self, ref, test, name):
        machine_tol = 1e-15

        if ref.shape != test.shape:
            self.error(f"    DIFF Different shape.\n  ref: {ref.shape}\n  test: {test.shape}")

        abs_ref = self.abs(ref)
        abs_test = self.abs(test)
        max_val = np.maximum(abs_ref, abs_test)
        not_small = max_val > machine_tol
        abs_diff = self.abs((ref - test)[not_small])
        max_val = max_val[not_small]
        abs_max_norm = np.max(abs_diff, initial=0)
        large_abs_diff = abs_diff > self.atol
        rel_diff = abs_diff[large_abs_diff] / max_val[large_abs_diff]
        rel_max_norm = np.max(rel_diff, initial=0)

        self.output.write(f"  {name}, adiff: {abs_max_norm}, rdiff: {rel_max_norm}")
        #if abs_max_norm > self.atol:
        #    n_abs_over = np.sum(abs_diff > self.atol)
        #    self.error(f"    DIFF {n_abs_over} values over atol: {self.atol}")
        if rel_max_norm > self.rtol:
            n_rel_over = np.sum(rel_diff > self.rtol)
            self.error(f"    DIFF {n_rel_over} values over rtol: {self.rtol} and atol: {self.atol}")

    @staticmethod
    def abs(array):
        if len(array.shape) == 1:
            return np.abs(array)
        elif len(array.shape) == 2:
            return np.linalg.norm(array, axis=1)
        else:
            return np.linalg.norm(array, axis=(1,2))


    def compare(self, ref_file, other_file, **kwargs):
        """
        Method can do anything as long as int value is returned
        :param reference_filepath:
        :param other_filepath:
        :param kwargs:
        :return: 1 - failed, 0 - success
        """
        self.atol = kwargs.get('atol', 1e-5)
        self.rtol = kwargs.get('rtol', 1e-2)
        self.interpolate = kwargs.get('interpolate', False)
        self.have_match = True

        #ref_file = "/home/jb/workspace/flow123d/tests/24_solute_dg/ref_out/01_sources/flow_test16/flow_test16-000000.vtu"
        #other_file = "/home/jb/workspace/flow123d/tests/24_solute_dg/test_results/01_sources.1/flow_test16/flow_test16-000000.vtu"
        ref_file = Paths.abspath(ref_file)
        other_file = Paths.abspath(other_file)
        
        self.compare_pipeline(ref_file, other_file)
        return not self.have_match




if __name__ == "__main__":
    vd = VTKDiff()
    a = get_parsed_args()
    if vd.compare(a.ref_file, a.test_file, atol=a.atol, rtol=a.rtol, interpolate=a.interpolate):
        # failed
        sys.exit(1)




"""
To use this checker, create check_rules section in config.yaml file
Files are read and per-line check is initiated.

Specify either regex keyword where regular expression is set
and/or
substr keyword which performs fast and simple sub string search.


On first match reading is interrupted.

e.g.:

common_config:
  proc: [1, 2]
  memory_limit: 1000
  check_rules:
    - regex:
        files: ["*.*"]
        regex: "[Ee]rror"           # matches the word error or Error in line
        substr: "error"             # looks for the word error in a line
"""


# shortcut to Regex so Python naming convention is not violated
vtkdiff = VTKDiff

