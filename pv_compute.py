from paraview.simple import *

# This is iscript for pvpython (Paraview Python) for computing
# L2 norm of difference of two output fields on different grids
#
# Usage: 
#
#   pvpython pv_compute.py  <coarse_file> <fine_file> <field_name>
#
# TODO: allow calculation of several fields in one run

import sys
import os 
import math

if len(sys.argv) < 3 :
  print "Missing argument! \n Usage: pv_compute.py <coarse_file> <fine_file> <field_name>\n"
  sys.exit(1)
  

coarse_input = sys.argv[1]
fine_input = sys.argv[2]
if not os.path.exists(coarse_input) :
  print "Missing file: ", coarse_input
  sys.exit(1)

if not os.path.exists(fine_input) :
  print "Missing file: ", fine_input
  sys.exit(1)
  
attribute=sys.argv[3]  


# create reader for legacy VTK files
# coarse datacoarse_input

coarse_sol = LegacyVTKReader(FileNames=coarse_input)
# fine data
fine_sol = LegacyVTKReader(FileNames=fine_input)

# Interpolate from coarse mesh to the fine mesh
# Input - values
# Source - mesh
coarse_on_fine = ResampleWithDataset( Source=fine_sol, Input=coarse_sol )

# compute difference of 'head' attribute (field)
diff = PythonCalculator( Input=[ coarse_on_fine, fine_sol ] )
diff.ArrayName = 'err_'+attribute;
dif_subexpr = "inputs[0].PointData['" + attribute + "'] - inputs[1].PointData['" + attribute + "']"
diff.Expression = "dot(" + dif_subexpr + "," + dif_subexpr + ")"

# compute square of the difference
#calc = Calculator(diff)
#calc.AttributeMode = 'point_data'
#calc.Function = diff.ArrayName + '*' + diff.ArrayName
#calc.ResultArrayName = 'err_' + attribute;

#Integrate square of the difference
iv = IntegrateVariables(diff)

# get output data from server
iv_output=paraview.servermanager.Fetch(iv)

# get Array with integrated value tuplets
array=iv_output.GetPointData().GetArray( diff.ArrayName )

# first value in the tuplet is what we need
print attribute + " L2 error: ", math.sqrt(array.GetTuple(0)[0])



