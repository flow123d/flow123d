"""
Script creates table of benchmark meshes
"""

import sys
import os
import gmsh
import pandas as pd


# initialize of environment and variables
directory = '../../benchmark_meshes'
gmsh.initialize()
table = pd.DataFrame( columns=['mesh_file', 'shape', 'uniformity', 'mesh_size', 'dimesnzion', 'number_of_elements'] )

# iterate over files in bench meshes directory
for file in os.listdir(directory):
    file_name, file_extension = os.path.splitext(file)
    
    if file_extension == '.msh':
        f = os.path.join(directory, file)
        mesh = gmsh.open(f)
        element_types, element_tags, node_tags = gmsh.model.mesh.getElements()
        n_elements = 0
        for i in range( len(element_tags) ):
            n_elements += len(element_tags[i])
        gmsh.clear()
        
        params = file_name.split("_")
        #square_2D_uniform_small
        table.loc[len(table.index)] = [file, params[0], params[2], params[3], params[1], n_elements]
    
print(table)