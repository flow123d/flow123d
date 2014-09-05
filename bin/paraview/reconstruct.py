#
# This script extract from given dataset all elements of specified type (lines, triangles, tetrahedra).
# The element type is set by the variable "select_type" below.
#
# TODO: pdo and pdi shares some (or all data), thus by with screw up input filter
# Make pdo NEW, pdo.Reset() is not enough.


import math
from ctypes import c_int

VTK_QUAD=9

TOL = 1e-5
X_SHIFT = 0.4

"""
Assumes a 2D mesh in XY plane with compatible 1D fracture on Y axis.
Duplicate nodes on fracture and shift elements in {x>0} by x_shift. 
Convert 1D line elements to rectangles x_shift wide. 
Convert cell and point data to the new mesh.
"""
def reconstruct_from_1d(pdi, pdo, x_shift):
    
    # Mark left=0 right=1, 1d elements=2
    # Mark their nodes.
    n_points_old=pdi.GetNumberOfPoints()
    n_cells_old=pdi.GetNumberOfCells()
    points_mask=[-1]*n_points_old
    cell_mask=[-1]*n_cells_old
    #cell_old_id=[]
    for i in range(n_cells_old):
        cell=pdi.GetCell(i)
        points=cell.GetPoints()
        x_center=0.0
        for i_point in range(points.GetNumberOfPoints()):
            point=points.GetPoint(i_point)
            x_center+=point[0]
        x_center/=points.GetNumberOfPoints()

        mask=-1
        if ( x_center < -TOL ):
            mask=0
        elif ( x_center > TOL ):
            mask=1
        else:
            mask=2
            
        cell_mask[i]=mask    
        for i_point in range(cell.GetNumberOfPoints()):
            id=cell.GetPointId(i_point)
            #print len(points_mask), "id:", id
            points_mask[id]=max( [points_mask[id], mask] )

    # Create new nodes (duplicate, shift) 
    # Make old_to_new_points map
    pdo.Reset()
    new_points = vtk.vtkPoints()
    point_new_id = [None]* n_points_old
    for i_old in range(n_points_old):
        if (points_mask[i_old]==0) :
            # just copy
            new_id=new_points.InsertNextPoint(pdi.GetPoint(i_old))
            point_new_id[i_old]=new_id
        elif (points_mask[i_old]==1) :
            # just shift
            point=pdi.GetPoint(i_old)
            point=(point[0]+x_shift, point[1], point[2])
            new_id=new_points.InsertNextPoint(point)
            point_new_id[i_old]=new_id            
        else :
            # duplicate
            point_b=point_a=pdi.GetPoint(i_old)
            point_b =(point_b[0]+x_shift, point_b[1], point_b[2])
            new_id_a=new_points.InsertNextPoint(point_a)
            new_id_b=new_points.InsertNextPoint(point_b)
            point_new_id[i_old]=[new_id_a, new_id_b]                        
    # Create new elements. (same numbering)
    # make cell list, substitute point id refs
    # Store Cell info to temporary structure, since it seems that pdi and pdo shares 
    # some storage.
    cell_info=[]
    for i_old in range(n_cells_old) :
        cell = pdi.GetCell(i_old)
        new_cell_point_ids = []
        if (cell_mask[i_old]==2) :
            # QUAD from LINE
            new_cell_point_ids.append(VTK_QUAD)
            points=[cell.GetPointId(0),cell.GetPointId(1)]
            new_cell_point_ids.append(point_new_id[ points[0] ][0])
            new_cell_point_ids.append(point_new_id[ points[1] ][0])
            new_cell_point_ids.append(point_new_id[ points[1] ][1])
            new_cell_point_ids.append(point_new_id[ points[0] ][1])
        else :
            # renumber points
            new_cell_point_ids.append(cell.GetCellType())
            for i_point in range(cell.GetNumberOfPoints()):
                id_old=cell.GetPointId(i_point)
                if (points_mask[id_old] == 2) :                   
                    cell_side = cell_mask[i_old]
                    new_point_id=point_new_id[ id_old ][cell_side]
                    new_cell_point_ids.append(new_point_id)
                else :
                    new_cell_point_ids.append( point_new_id[ id_old ] )
        cell_info.append(new_cell_point_ids)
    # Store new mesh into pdo
    pdo.SetPoints(new_points)
    for cell_list in cell_info :
        cell_type=cell_list[0]
        id_list=vtk.vtkIdList()
        for point_id in cell_list[1:] :
            id_list.InsertNextId(point_id)
        pdo.InsertNextCell( cell_type, id_list)
    # Reuse CellData
    for i in range(pdi.GetCellData().GetNumberOfArrays()) :
        print pdi.GetCellData().GetArray(i).GetName()
        pdo.GetCellData().AddArray( pdi.GetCellData().GetArray(i) )           
    # Convert Point data, For every point put its value 
    # to new location.
    # make point data arrays
    for i in range(pdi.GetPointData().GetNumberOfArrays()) :
        oldArray=pdi.GetPointData().GetArray(i)
        print oldArray.GetName()
        newArray = vtk.vtkDoubleArray()
        newArray.SetName(oldArray.GetName())
        newArray.SetNumberOfComponents(oldArray.GetNumberOfComponents())
        newArray.SetNumberOfTuples( new_points.GetNumberOfPoints() )
        
        
        for i_old in range( n_points_old ) :
            if (points_mask[i_old]==2) :
                i_new = point_new_id[ i_old ]  
                newArray.SetTupleValue(i_new[0], oldArray.GetTuple(i_old) )
                newArray.SetTupleValue(i_new[1], oldArray.GetTuple(i_old) )
            else:    
                i_new = point_new_id[ i_old ]  
                newArray.SetTupleValue(i_new, oldArray.GetTuple(i_old) )        
        pdo.GetPointData().AddArray(newArray)

################################################
# The filter with appropriate predicate.

pdi = self.GetPolyDataInput()
pdo =self.GetOutput()
pdo.DeepCopy(pdi)
reconstruct_from_1d(pdi, pdo, X_SHIFT)

