"""
Script to extract given set of elements from the given mesh file.
Renumber nodes starting from 1 and numbering nodes of elements according to order given at input.
Usage:

    python extract_elements.py <mesh_file> <element id list>
"""

import sys

ele_sizes = {
    23 : 0, # point
    1 : 2,  # line
    2 : 3,  # triangle
    4 : 4  # tetrahedron
}

msh_file = sys.argv[1]
ele_ids=[]
for ie in sys.argv[2:]:
    ele_ids.append( int(ie) )
    
header=[]    
elements=[]
nodes=[]
nodes_in_use={}
old_elements={}
n_nodes=0



with open(msh_file, "r") as mf:
    for line in mf:
        if '$Nodes' in line:
            break
        header.append(line)

    for line in mf:
        if '$Elements' in line:
            break
        
    # extract elements
    n_elements=int(mf.readline())
    for i_ele in range(n_elements):
        line = mf.readline().split()
        if int(line[0]) in ele_ids:
           ele_items = [ int(i) for i in line ]
           old_elements[ ele_items[0] ] = ele_items
    
    # convert elements
    for ele_id in ele_ids:
        ele_items = old_elements[ele_id]
        n_vtx = ele_items[1]
        n_tags = ele_items[2]
        ele_nodes = ele_items[3+n_tags : ]
        # collect and rename nodes
        ele_new_nodes =[]
        for node in ele_nodes:
            if not node in nodes_in_use:
                n_nodes+=1
                nodes_in_use[node] = n_nodes
            new_node = nodes_in_use[node]   
            ele_new_nodes.append(new_node)
        # convert ID and nodes    
        ele_items[0] = len(elements) + 1
        ele_items[3+n_tags : ] = ele_new_nodes    
        ele_str = " ".join([str(x) for x in ele_items])
        elements.append( ele_str )          
    
    mf.seek(0)

    # extract and convert nodes
    for line in mf:
        if '$Nodes' in line:
            break
    
    
    # extract and convert elements
    n_nodes=int(mf.readline())
    for i_node in range(n_nodes):
        line = mf.readline().split()
        if int(line[0]) in nodes_in_use.keys():
           line[0] = str( nodes_in_use[ int(line[0]) ] )        
           nodes.append(" ".join(line))
           

def write_block(omf, name, lines):
    print("$"+name, file=omf)
    print("%d"%(len(lines)), file=omf)
    for l in lines:
        print(l, file=omf)
    print("$End" + name, file=omf)
    
           
with open("sub_"+msh_file, "w") as omf:
    for l in header:
        omf.write(l)
    write_block(omf, "Nodes", nodes)    
    write_block(omf, "Elements", elements)
    