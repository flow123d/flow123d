import itertools

# 1D: given by vertices
t1d = [0, 1]

# 2D: given by sides
t2d = [[0,1], [0,2], [1,2]]

# 3D: given by faces
t3d = [[0,1,2], [0,1,3], [0,2,3], [1,2,3]]

edges = [[0,1], [0,2], [0,3], [1, 2], [1,3], [2,3]]

def find_edge(xedge):
    for ie, e in enumerate(edges):
        if xedge == e:
            return (ie, True)
        elif xedge == [e[1], e[0]]:
            return (ie, False)
    assert(False)

def face_edges(face):
     f_edges = [[face[s0], face[s1]] for s0, s1 in t2d]
     return [find_edge(e)[0] for e in f_edges]
  
t3d_face_edges = [face_edges(face) for face in t3d]

# 64 possible edge orientations of a 3d tetrahedron
# try to find vertex permutation so that it match the native orientation

def bitfield(n):
     # bit field for 64 bit integer
     return [n >> i & 1 for i in range(5,-1,-1)]

def bitnumber(field):
    n = 0
    for f in field:
        n += f
        n = n << 1
    return n >> 1
  
assert bitnumber(bitfield(42)) == 42
assert bitnumber(bitfield(13)) == 13

def permute(x, perm):
    return [x[i] for i in perm]

assert ['B', 'C', 'A'] == permute(['A', 'B', 'C'], (1,2,0))
  

def is_sorted(x):
    return int(x[0] < x[1])
             
def side_orientations(perm_2d):
    return [is_sorted(permute(side, perm_2d)) for side in t2d]


valid_face_edge_ori = [1,1,0,1,1,0,1,1]
  
def is_possible(edge_ori):
    # first 3 faces must be valid
    for i_face in range(3):
        face_edge_ori = [edge_ori[ie] for ie in t3d_face_edges[i_face]]
        if not valid_face_edge_ori[bitnumber(face_edge_ori)]:
            return False
    return True

def edge_orientations(elem_perm):
    vertices = permute(list(range(4)), elem_perm)
    ori = [int(vertices[v0] < vertices[v1]) for v0, v1 in edges]
    #print("   ", ori)
    return ori

def valid_edge_orientations_on_element():
    for i_edge_orientations in range(64):
        edge_ori = bitfield(i_edge_orientations)
        if is_possible(edge_ori):
            valid_el_perms = [elem_perm for elem_perm in itertools.permutations(list(range(4))) if edge_ori == edge_orientations(elem_perm)]
            msg = str(valid_el_perms)
            #msg = "possible"
        else:
            msg = "impossible"
        print(edge_ori, msg)
            
            
def node_permutations():
    """
    Generate a C++ code for the array providing the new_to_old mapping for all 64 bitfield values of the 
    all six node comparisons (assumes 4 nodes for all elements).
    """
    comparisons = 64 * [0]
    # Generate permutations as new to old node mapping for single element.
    for elem_perm in itertools.permutations(list(range(4))):
        old_to_new = list(range(4))
        for new, old in enumerate(elem_perm):
            old_to_new[old] = new
        perm_index = 0
        for i, edge in enumerate(edges):
            n0, n1 = edge
            edge_comparison = (old_to_new[n0] > old_to_new[n1])
            perm_index *= 2
            perm_index += edge_comparison
    
        out_list = ",".join(map(str, elem_perm))    
        print(f"element_nodes_original[{perm_index}] = {{{out_list}}};")
            
        
node_permutations()
#list()
