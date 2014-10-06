#
# This script extract from given dataset all elements of specified type (lines, triangles, tetrahedra).
# The element type is set by the variable "select_type" below.
#

import math
from ctypes import c_int

# select_type
#
# Choose VTK elements of given type. According to VTK file format documentation
# (http://www.vtk.org/VTK/img/file-formats.pdf), the following values are possible:
# 
# 1 VTK_VERTEX
# 3 VTK_LINE
# 5 VTK_TRIANGLE
# 10 VTK_TETRA



# Predicates for cell selection. Resulting mesh will contain only cells
# for which predicate returns true.
def is_line(cell):
    return cell.GetCellType() == 3
  
def is_triangle(cell):
    return cell.GetCellType() == 5

def is_XY_plane(cell):
    points=cell.GetPoints()
    z_center=0.0
    for i_point in range(points.GetNumberOfPoints()):
        point=points.GetPoint(i_point)
        z_center+=point[2]
    z_center/=points.GetNumberOfPoints()
    return abs(z_center) < 1e-5


def submesh(pdi, predicate, pdo):
    # select new cells, cell_old_id
    # mark necessary points,
    n_points_old=pdi.GetNumberOfPoints()
    points_mask=[0]*n_points_old;
    cell_old_id=[]
    for i in range(pdi.GetNumberOfCells()):
        cell=pdi.GetCell(i)
        if ( predicate(cell) ):
            cell_old_id.append(i)
            for i_point in range(cell.GetNumberOfPoints()):
                points_mask[cell.GetPointId(i_point)]=1

    pdo.Reset()
    # make point list, point_new_id
    new_points = vtk.vtkPoints()
    point_new_id = [-1]* n_points_old
    for i_old in range(n_points_old):
        if points_mask[i_old] :
            new_id=new_points.InsertNextPoint(pdi.GetPoint(i_old))
            point_new_id[i_old]=new_id
    pdo.SetPoints(new_points)
    
    # make cell list, substitute point id refs
    for i_old in cell_old_id :
        cell = pdi.GetCell(i_old)
        new_cell_point_ids = vtk.vtkIdList()
        for i_point in range(cell.GetNumberOfPoints()):
            new_cell_point_ids.InsertNextId( point_new_id[ cell.GetPointId(i_point) ] )
        pdo.InsertNextCell( cell.GetCellType(), new_cell_point_ids); # InsertNextCell( int  type, vtkIdList *  ptIds );"""

    # make cell data arrays
    for i in range(pdi.GetCellData().GetNumberOfArrays()) :
        oldArray=pdi.GetCellData().GetArray(i)
        newArray = vtk.vtkDoubleArray()
        newArray.SetName(oldArray.GetName())
        newArray.SetNumberOfComponents(oldArray.GetNumberOfComponents())
        newArray.SetNumberOfTuples(cell_old_id.__len__())
        
        i_new=0
        for i_old in cell_old_id :
            newArray.SetTupleValue(i_new, oldArray.GetTuple(i_old) )
            i_new+=1
            
        pdo.GetCellData().AddArray(newArray)

    # make point data arrays
    for i in range(pdi.GetPointData().GetNumberOfArrays()) :
        oldArray=pdi.GetPointData().GetArray(i)
        newArray = vtk.vtkDoubleArray()
        newArray.SetName(oldArray.GetName())
        newArray.SetNumberOfComponents(oldArray.GetNumberOfComponents())
        newArray.SetNumberOfTuples( new_points.GetNumberOfPoints() )
        
        
        for i_old in range( len( point_new_id ) ) :
            i_new = point_new_id[ i_old ]
            if ( i_new != -1 ) :
                newArray.SetTupleValue(i_new, oldArray.GetTuple(i_old) )
        
        pdo.GetPointData().AddArray(newArray)

################################################
# The filter with appropriate predicate.

pdi = self.GetPolyDataInput()
pdo = self.GetOutput()
submesh(pdi, is_XY_plane, pdo)
