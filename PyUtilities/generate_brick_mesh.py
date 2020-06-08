import numpy as np
import h5py

if __name__ == "__main__":
    h5_file = h5py.File("..\\Asset\\brick_mesh.h5", "w")
    mesh_grp = h5_file.create_group("Mesh")
    # write nodes
    nodes = [
        [2, -1.0, -1.0, -1.0],
        [3,  1.0, -1.0, -1.0],
        [4,  1.0,  1.0, -1.0],
        [5, -1.0,  1.0, -1.0],
        [6, -1.0, -1.0,  1.0],
        [7,  1.0, -1.0,  1.0],
        [8,  1.0,  1.0,  1.0],
        [9, -1.0,  1.0,  1.0]
    ]
    node_type = np.dtype([("index", np.int64), ("x", np.float64), ("y", np.float64), ("z", np.float64)])
    nd = np.empty(len(nodes), dtype=node_type)
    for i in range(len(nodes)):
        nd[i]["index"] = nodes[i][0]
        nd[i]["x"] = nodes[i][1]
        nd[i]["y"] = nodes[i][2]
        nd[i]["z"] = nodes[i][3]
    n_dset = mesh_grp.create_dataset("Nodes", data=nd)
    # write elements
    elems = [
        [3, 2, 3, 5, 6],
        [4, 3, 4, 5, 8],
        [5, 3, 5, 7, 6],
        [7, 3, 5, 7, 8],
        [8, 6, 7, 9, 5],
        [1, 7, 8, 9, 5]
    ]
    elem_type = np.dtype([("index", np.int64), ("n1", np.int64), ("n2", np.int64), ("n3", np.int64), ("n4", np.int64)])
    ed = np.empty(len(elems), dtype=elem_type)
    for i in range(len(elems)):
        ed[i]["index"] = elems[i][0]
        ed[i]["n1"] = elems[i][1]
        ed[i]["n2"] = elems[i][2]
        ed[i]["n3"] = elems[i][3]
        ed[i]["n4"] = elems[i][4]
    e_dset = mesh_grp.create_dataset("Elements", data=ed)
    h5_file.close()
