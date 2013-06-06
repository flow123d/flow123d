#
# This script extract from given dataset all elements of specified type (lines, triangles, tetrahedra).
# The element type is set by the variable "select_type" below.
#

import math

# select_type
#
# Choose VTK elements of given type. According to VTK file format documentation
# (http://www.vtk.org/VTK/img/file-formats.pdf), the following values are possible:
# 
# 1 VTK_VERTEX
# 3 VTK_LINE
# 5 VTK_TRIANGLE
# 10 VTK_TETRA

select_type=5

pdi = self.GetPolyDataInput()
pdo = self.GetOutput()


# make filter list
new_to_old=[]
for i in range (pdi.GetNumberOfCells()):
      if pdi.GetCellType(i)== select_type:
          new_to_old.append(i)   
          
# make new cell list
newcells = vtk.vtkCellArray()
for i_old in new_to_old :
    cell = pdi.GetCell(i_old)
    newcells.InsertNextCell(cell)

pdo.SetCells( select_type, newcells )

# make cell data arrays
for i in range(pdi.GetCellData().GetNumberOfArrays()) :
    oldArray=pdi.GetCellData().GetArray(i)
    newArray = vtk.vtkDoubleArray()
    newArray.SetName(oldArray.GetName())
    newArray.SetNumberOfComponents(oldArray.GetNumberOfComponents())
    newArray.SetNumberOfTuples(new_to_old.__len__())
    
    i_new=0
    for i_old in new_to_old :
        newArray.SetTupleValue(i_new, oldArray.GetTuple(i_old) )
        i_new+=1
        
    pdo.GetCellData().AddArray(newArray)

# can not simply copy Point Data, that needs modification of points numbering    
