import math

select_type=3

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
