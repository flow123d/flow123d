# Usage: python3 compare_view.py FILE
# 
# For given file path finds both reference and the test file and starts two instances of corresponding vizualization tool GMSH or Paraview.
# file can be either the reference path, e.g. 02_dirichlet/transport/transport-000001.vtu
# or test path with number of processors, e.g. 02_dirichlet.2/transport/transport-000001.vtu
#
# After one viewer is closed the script automaticaly kills the other one.


import subprocess
import argparse
import re
import os
import concurrent.futures
import time


gmsh = "/home/jb/bin/gmsh4"
paraview = "/home/jb/local/ParaView-5.4.0-Qt5-OpenGL2-MPI-Linux-64bit/bin/paraview"


class ProcessThread:
    class NoneResult:
        pass
    none_result = NoneResult()
    
    executor = concurrent.futures.ThreadPoolExecutor(max_workers=2)
    
    def __init__(self, command, **kwargs):
        """
        command  - list containting the command ant its args as passed to the Popen
        """
        self.command = command
        self.future = ProcessThread.executor.submit(self.cmd_call, command, **kwargs)

    def cmd_call(self, command, **kwargs):
        self.sub_process = subprocess.Popen(command, **kwargs)
        self.sub_process.wait()
        #print(f"END CMD CALL: {command}")
        return 
    
    def result(self):
        try:
            return self.future.result(timeout=0)
        except concurrent.futures.TimeoutError:
            return ProcessThread.none_result
    
    def finished(self):
        res = self.result()
        #print(f"{res} ... for {self.command}")
        return res != ProcessThread.none_result
        
          
    def kill(self):
        """
        Kill the subprocess.
        """
        if not self.future.cancel():
            if self.future.running():
                self.sub_process.kill()
            else:
                return
            
    


def get_parsed_args():
    parser = argparse.ArgumentParser(
        prog="vtk_diff",
        description="Numerical comparison of two VTK/VTU files: reference vs. test result.")
    parser.add_argument('ref_path',
                        help="Reference result path relative to ref_out.")

    a = parser.parse_args()
    print(a)
    return a
  
def path_iter(path):
    while path:
      path, tail = os.path.split(path)
      yield tail
    
    
def path_split(path):
    rev =  [f for f in path_iter(path)]
    rev.reverse()
    return rev


def find_test_path(ref_path):
    ref_root = 'ref_out'
    test_root = 'test_results'
    path_list = path_split(ref_path)
    
    base, ext = os.path.splitext(path_list[0])
    if not ext:
        ext = ".1"        
    ref_file = os.path.join(ref_root, base, *path_list[1:])
    test_file = os.path.join(test_root, base+ext, *path_list[1:])
    assert os.path.isfile(ref_file), f"not found: {ref_file}"
    assert os.path.isfile(test_file), f"not found: {test_file}"
    return ref_file, test_file
       
    
  
  
def view(path):
    _, ext = os.path.splitext(path)
    if ext == ".msh":
        return ProcessThread([gmsh, path])
    else:
        return ProcessThread([paraview, path])

    
if __name__ == "__main__":    
    a = get_parsed_args()
    
    ref_path,  test_path = find_test_path(a.ref_path)
    ref_view = view(ref_path)
    test_view = view(test_path)
    while not ref_view.finished() and not test_view.finished():
        time.sleep(1)
    ref_view.kill()
    test_view.kill()
            
            
    
    
