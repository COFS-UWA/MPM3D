import math
import numpy as np

def generate_cone_surf_mesh(num, filename):
    node_coord = [[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]]
    elem_index = []
    dtheta = math.pi * 2.0 / float(num)
    for i in range(num):
        node_coord.append([math.cos(float(i) * dtheta), math.sin(float(i) * dtheta), 1.0])
    for i in range(num):
        node_coord.append([math.cos(float(i) * dtheta), math.sin(float(i) * dtheta), 0.0])
    node_num = 2 + 2 * num
    for i in range(num-1):
        elem_index.append([0, i+2, i+3])
    elem_index.append([0, num+1, 2])
    for i in range(num-1):
        elem_index.append([1, num+i+3, num+i+2])
    elem_index.append([1, num+2, 2*num+1])
    for i in range(num-1):
        elem_index.append([i+2, num+i+2, num+i+3])
        elem_index.append([i+2, num+i+3, i+3])
    elem_index.append([num+1, 2*num+1, num+2])
    elem_index.append([num+1, num+2, 2])
    elem_num = 4 * num
    
    with open(filename, "w") as cfile:
        cfile.write("""static const unsigned int cone_node_num = %d;
static const unsigned int cone_elem_num = %d;

static const float cone_nodes[] = {
""" % (node_num, elem_num))
        
        for n_id in range(node_num-1):
            node = node_coord[n_id]
            cfile.write("    %ff, %ff, %ff, // node %d\n" % (node[0], node[1], node[2], n_id))
        node = node_coord[node_num-1]
        cfile.write("    %ff, %ff, %ff // node %d\n};" % (node[0], node[1], node[2], node_num-1))
        
        cfile.write("""
        
static const unsigned int cone_elems[] = {
""")
        
        for e_id in range(elem_num-1):
            elem = elem_index[e_id]
            cfile.write("    %d, %d, %d, // elem %d\n" % (elem[0], elem[1], elem[2], e_id))
        elem = elem_index[elem_num-1]
        cfile.write("    %d, %d, %d // elem %d\n};" % (elem[0], elem[1], elem[2], elem_num-1))

if __name__ == "__main__":
    generate_cone_surf_mesh(120, "cone_surf_mesh.h")
