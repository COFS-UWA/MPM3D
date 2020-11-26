import math
import numpy as np

def generate_cylinder_surf_mesh(num, filename):
    node_coord = [[0.0, 0.0, 0.5], [0.0, 0.0, -0.5]]
    elem_index = []
    dtheta = math.pi * 2.0 / float(num)
    for i in range(num):
        node_coord.append([math.cos(float(i) * dtheta), math.sin(float(i) * dtheta), 0.5])
    for i in range(num):
        node_coord.append([math.cos(float(i) * dtheta), math.sin(float(i) * dtheta), -0.5])
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
        cfile.write("""#ifndef __Cylinder_Surf_Mesh_Data__
#define __Cylinder_Surf_Mesh_Data__

static const unsigned int cylinder_node_num = %d;

static const float cylinder_nodes[] = {
""" % node_num)
        
        for n_id in range(node_num-1):
            node = node_coord[n_id]
            cfile.write("    %ff, %ff, %ff, // node %d\n" % (node[0], node[1], node[2], n_id))
        node = node_coord[node_num-1]
        cfile.write("    %ff, %ff, %ff // node %d\n};" % (node[0], node[1], node[2], node_num-1))
        
        cfile.write("""

static const unsigned int cylinder_elem_num = %d;

static const unsigned int cylinder_elems[] = {
""" % elem_num)
        
        for e_id in range(elem_num-1):
            elem = elem_index[e_id]
            cfile.write("    %d, %d, %d, // elem %d\n" % (elem[0], elem[1], elem[2], e_id))
        elem = elem_index[elem_num-1]
        cfile.write("    %d, %d, %d // elem %d\n};" % (elem[0], elem[1], elem[2], elem_num-1))

        cfile.write("""

#endif""")

if __name__ == "__main__":
    generate_cylinder_surf_mesh(120, "cylinder.h")
