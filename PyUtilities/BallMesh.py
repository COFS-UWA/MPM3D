import math
import numpy as np

def move_to_sphere(vec3):
    vec_len = math.sqrt(vec3[0] * vec3[0] + vec3[1] * vec3[1] + vec3[2] * vec3[2])
    vec3[0] /= vec_len
    vec3[1] /= vec_len
    vec3[2] /= vec_len

# Tessellate sphere into triangle mesh
class BallMesh:
    def __init__(self):
        self.node_num = 0
        self.nodes = None
        self.elem_num = 0
        self.elems = None
    
    def gen_regular_icosahedron(self):
        # nodes
        m = math.sqrt(50.0 - 10.0 * math.sqrt(5.0)) / 10.0
        n = math.sqrt(50.0 + 10.0 * math.sqrt(5.0)) / 10.0
        self.node_num = 12
        self.nodes = np.array([
            [ m, 0.0,  n],
            [ m, 0.0, -n],
            [-m, 0.0,  n],
            [-m, 0.0, -n],
            [0.0,  n,  m],
            [0.0, -n,  m],
            [0.0,  n, -m],
            [0.0, -n, -m],
            [ n,  m, 0.0],
            [-n,  m, 0.0],
            [ n, -m, 0.0],
            [-n, -m, 0.0]
        ])
        # elems
        self.elem_num = 20
        self.elems = np.array([
            [6, 4, 8],
            [9, 4, 6],
            [6, 3, 9],
            [6, 1, 3],
            [6, 8, 1],
            [8, 10, 1],
            [8, 0, 10],
            [8, 4, 0],
            [4, 2, 0],
            [4, 9, 2],
            [9, 11, 2],
            [9, 3, 11],
            [3, 1, 7],
            [1, 10, 7],
            [10, 0, 5],
            [0, 2, 5],
            [2, 11, 5],
            [3, 7, 11],
            [5, 11, 7],
            [10, 5, 7]
        ], dtype=np.uint32)
    
    def subdivide(self, other):
        self.node_num = other.node_num + other.elem_num * 3
        self.nodes = np.zeros([self.node_num, 3])
        self.elem_num = other.elem_num * 4
        self.elems = np.zeros([self.elem_num, 3], dtype=np.uint32)
        # copy data of old node
        for n_id in range(other.node_num):
            o_node = other.nodes[n_id]
            n_node = self.nodes[n_id]
            n_node[0] = o_node[0]
            n_node[1] = o_node[1]
            n_node[2] = o_node[2]
        for e_id in range(other.elem_num):
            o_elem = other.elems[e_id]
            n_elem = self.elems[e_id]
            n_elem[0] = o_elem[0]
            n_elem[1] = o_elem[1]
            n_elem[2] = o_elem[2]
        # generate data for new node
        nid_off = other.node_num
        eid_off = other.elem_num
        for i in range(other.elem_num):
            elem = self.elems[i]
            node1 = self.nodes[elem[0]]
            node2 = self.nodes[elem[1]]
            node3 = self.nodes[elem[2]]
            n_node1 = self.nodes[nid_off + 3*i]
            n_node1[0] = (node1[0] + node2[0]) * 0.5
            n_node1[1] = (node1[1] + node2[1]) * 0.5
            n_node1[2] = (node1[2] + node2[2]) * 0.5
            move_to_sphere(n_node1)
            n_node2 = self.nodes[nid_off + 3*i+1]
            n_node2[0] = (node2[0] + node3[0]) * 0.5
            n_node2[1] = (node2[1] + node3[1]) * 0.5
            n_node2[2] = (node2[2] + node3[2]) * 0.5
            move_to_sphere(n_node2)
            n_node3 = self.nodes[nid_off + 3*i+2]
            n_node3[0] = (node3[0] + node1[0]) * 0.5
            n_node3[1] = (node3[1] + node1[1]) * 0.5
            n_node3[2] = (node3[2] + node1[2]) * 0.5
            move_to_sphere(n_node3)
            n_elem1 = self.elems[eid_off + 3*i]
            n_elem1[0] = elem[0]
            n_elem1[1] = nid_off + 3*i
            n_elem1[2] = nid_off + 3*i+2
            n_elem2 = self.elems[eid_off + 3*i+1]
            n_elem2[0] = elem[1]
            n_elem2[1] = nid_off + 3*i+1
            n_elem2[2] = nid_off + 3*i
            n_elem3 = self.elems[eid_off + 3*i+2]
            n_elem3[0] = elem[2]
            n_elem3[1] = nid_off + 3*i+2
            n_elem3[2] = nid_off + 3*i+1
            elem[0] = nid_off + 3*i
            elem[1] = nid_off + 3*i+1
            elem[2] = nid_off + 3*i+2
    
    def write_to_cfile(self, filename):
        with open(filename, "w") as cfile:
            cfile.write("""#ifndef __Ball_Mesh_Data__
#define __Ball_Mesh_Data__

static const unsigned int ball_node_num = %d;

static const float ball_nodes[] = {
""" % self.node_num)
            
            for n_id in range(self.node_num-1):
                node = self.nodes[n_id]
                cfile.write("    %f, %f, %f, // node %d\n" % (node[0], node[1], node[2], n_id))
            node = self.nodes[self.node_num-1]
            cfile.write("    %f, %f, %f // node %d\n};" % (node[0], node[1], node[2], self.node_num-1))
            
            cfile.write("""

static const unsigned int ball_elem_num = %d;

static const unsigned int ball_elems[] = {
""" % self.elem_num)

            for e_id in range(self.elem_num-1):
                elem = self.elems[e_id]
                cfile.write("    %d, %d, %d, // elem %d\n" % (elem[0], elem[1], elem[2], e_id))
            elem = self.elems[self.elem_num-1]
            cfile.write("    %d, %d, %d // elem %d\n};" % (elem[0], elem[1], elem[2], e_id))

            cfile.write("""

#endif""")

if __name__ == "__main__":
    ball0 = BallMesh()
    ball0.gen_regular_icosahedron()
    ball1 = BallMesh()
    ball1.subdivide(ball0)
    ball1.write_to_cfile("ball_mesh_data.h")
    #ball2 = BallMesh()
    #ball2.subdivide(ball) 