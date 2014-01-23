from paraview.simple import *

# This is script for pvpython (Paraview Python) for computing
# L2 norm of difference of two output fields on possibly different grids
#
# Usage: 
#
#   pvpython pv_compute.py  <coarse_file> <fine_file> <field_name1> <field_name2> ...
#
# TODO: allow calculation of several fields in one run
#       check type of data and produce same in python calculator, check existence of field names
#
# coarse_file - either VTK file or python script (extension.py) to compute element data fields on fine_file grid 
# fine_file - VTK file with reference grid and data
# field_names - list of fields to differs

import sys
import os 
import math
import re

#"""
if len(sys.argv) < 3 :
  print "Missing argument! \n Usage: pv_compute.py <coarse_file> <fine_file> <field_name>\n"
  sys.exit(1)
  

coarse_input = sys.argv[1]
fine_input = sys.argv[2]
if not os.path.exists(coarse_input) :
  print "Missing file: ", coarse_input
  sys.exit(1)

if re.search("\.py$",coarse_input) :
  script_on_coarse=True
else :
  script_on_coarse=False
  
if not os.path.exists(fine_input) :
  print "Missing file: ", fine_input
  sys.exit(1)
  
attribute_list=sys.argv[3:]  

#"""
"""
coarse_input="/home/jb/Projekty/09_02_pukliny_teorie/11_04_numericke testy/noncompatible_2/analytical_solution.py"
script_on_coarse=True
fine_input="/home/jb/Projekty/09_02_pukliny_teorie/11_04_numericke testy/noncompatible_2/dx_3_r_1_s_100/square.vtu"
attribute_list=["element_scalars", "element_vectors"]
"""


#fine_sol = LegacyVTKReader(FileNames=fine_input)
fine_sol = XMLUnstructuredGridReader(FileName=fine_input)

if not script_on_coarse :
  # coarse_sol = LegacyVTKReader(FileNames=coarse_input)
  coarse_sol = XMLUnstructuredGridReader(FileName=coarse_input)

  # Interpolate from coarse mesh to the fine mesh
  # Input - values
  # Source - mesh
  coarse_on_fine = ResampleWithDataset( Source=fine_sol, Input=coarse_sol )
else :
  
  # create Programmable Filter with script
  script_file=open(coarse_input, "r")
  script_str=script_file.read()
  script_file.close()
  coarse_on_fine = ProgrammableFilter( Input=fine_sol, CopyArrays=0, Script=script_str)

print "Sources created.\n"
  
  
py_calc_list=[]
for attribute in attribute_list :
  # compute difference of 'head' attribute (field)
  diff = PythonCalculator( Input=[ coarse_on_fine, fine_sol ] )
  diff.ArrayName = 'err_'+attribute;
  dif_subexpr = "inputs[0].CellData['" + attribute + "'] - inputs[1].CellData['" + attribute + "']"
  diff.Expression = "dot(" + dif_subexpr + "," + dif_subexpr + ")"
  diff.ArrayAssociation=1 	# Cell data
  diff.CopyArrays=0		# we want only err attributes in 'all_err_attributes'
  py_calc_list.append(diff)
  
all_err_attributes= AppendAttributes(Input=py_calc_list)


#Integrate square of the difference
iv_2d = IntegrateVariables(all_err_attributes)

# Filter 1D data and integrate on them
script_file=open("../submesh_1d.py", "r")
script_str=script_file.read()
script_file.close()
pf_1d = ProgrammableFilter( Input=all_err_attributes, CopyArrays=0, Script=script_str)
iv_1d = IntegrateVariables(pf_1d)

# stop the program and provide console; continue Ctrl-D
#import code; code.interact(local=locals())

print "Updating pipeline..."
iv_2d.UpdatePipeline()
iv_1d.UpdatePipeline()
print "done. \n Fetching Integrated data.."

# get output data from server
iv_2d_output=paraview.servermanager.Fetch(iv_2d)
iv_1d_output=paraview.servermanager.Fetch(iv_1d)

for attribute in attribute_list :
  # get Array with integrated value tuplets
  array=iv_2d_output.GetCellData().GetArray( 'err_'+attribute )
  # first value in the tuplet is what we need
  print attribute + " 2D L2 error: ", math.sqrt(array.GetTuple(0)[0])

for attribute in attribute_list :
  # get Array with integrated value tuplets
  array=iv_1d_output.GetCellData().GetArray( 'err_'+attribute )
  # first value in the tuplet is what we need
  print attribute + " 1D L2 error: ", math.sqrt(array.GetTuple(0)[0])

print "Finished."  


