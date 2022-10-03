import numpy as np
import h5py

def generate_rect_mesh_1byn(grid_num, grid_size):
    # generate nodes
    nodes = [
        [0, 0.0, 0.0],
        [1, 0.0, grid_size]
    ]
    for i in range(grid_num):
        x = float(1 + i) * grid_size
        nodes.append([2*i+2, x, 0.0])
        nodes.append([2*i+3, x, grid_size])
    # generate elems
    elems = []
    for i in range(grid_num):
        elems.append([2*i, i*2, i*2+2, i*2+3])
        elems.append([2*i+1, i*2, i*2+3, i*2+1])
    return nodes, elems
    
if __name__ == "__main__":
    grid_size = 0.06
    grid_num = 5
    nodes, elems = generate_rect_mesh_1byn(grid_num, grid_size)
    h5_file = h5py.File("../Asset/rect_mesh_1by%d.h5" % grid_num, "w")
    mesh_grp = h5_file.create_group("Mesh")
    # write nodes
    node_type = np.dtype([("index", np.int64), ("x", np.float64), ("y", np.float64)])
    node_num = len(nodes)
    nd = np.empty(node_num, dtype = node_type)
    for i in range(node_num):
        nd[i]["index"] = nodes[i][0]
        nd[i]["x"] = nodes[i][1]
        nd[i]["y"] = nodes[i][2]
    n_dset = mesh_grp.create_dataset("Nodes", data = nd)
    # write elements
    elem_type = np.dtype([("index", np.int64), ("n1", np.int64), ("n2", np.int64), ("n3", np.int64)])
    elem_num = len(elems)
    ed = np.empty(elem_num, dtype = elem_type)
    for i in range(elem_num):
        ed[i]["index"] = elems[i][0]
        ed[i]["n1"] = elems[i][1]
        ed[i]["n2"] = elems[i][2]
        ed[i]["n3"] = elems[i][3]
    e_dset = mesh_grp.create_dataset("Elements", data = ed)
    h5_file.close()
