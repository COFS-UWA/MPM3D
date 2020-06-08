import numpy as np
import h5py

if __name__ == "__main__":
    level_num = 5
    elem_size = 0.1
    h5_file = h5py.File("..\\Asset\\bar_mesh1.h5", "w")
    mesh_grp = h5_file.create_group("Mesh")
    # write nodes
    nodes = []
    for i in range(level_num+1):
        i4 = 4 * i
        h = float(i) * elem_size
        nodes.append((0+i4, 0.0,       0.0,       h))
        nodes.append((1+i4, elem_size, 0.0,       h))
        nodes.append((2+i4, elem_size, elem_size, h))
        nodes.append((3+i4, 0.0,       elem_size, h))
    node_type = np.dtype([
        ("index", np.int64),
        ("x", np.float64),
        ("y", np.float64),
        ("z", np.float64)
        ])
    nd = np.empty(len(nodes), dtype=node_type)
    for i in range(len(nodes)):
        nd[i]["index"] = nodes[i][0]
        nd[i]["x"] = nodes[i][1]
        nd[i]["y"] = nodes[i][2]
        nd[i]["z"] = nodes[i][3]
    n_dset = mesh_grp.create_dataset("Nodes", data=nd)
    # write elements
    elems = []
    for i in range(level_num):
        i4 = 4 * i
        i5 = 5 * i
        elems.append((0+i5, 0+i4, 1+i4, 3+i4, 4+i4))
        elems.append((1+i5, 1+i4, 2+i4, 3+i4, 6+i4))
        elems.append((2+i5, 4+i4, 7+i4, 6+i4, 3+i4))
        elems.append((3+i5, 4+i4, 6+i4, 5+i4, 1+i4))
        elems.append((4+i5, 1+i4, 4+i4, 6+i4, 3+i4))
    elem_type = np.dtype([
        ("index", np.int64),
        ("n1", np.int64),
        ("n2", np.int64),
        ("n3", np.int64),
        ("n4", np.int64)
        ])
    ed = np.empty(len(elems), dtype=elem_type)
    for i in range(len(elems)):
        ed[i]["index"] = elems[i][0]
        ed[i]["n1"] = elems[i][1]
        ed[i]["n2"] = elems[i][2]
        ed[i]["n3"] = elems[i][3]
        ed[i]["n4"] = elems[i][4]
    e_dset = mesh_grp.create_dataset("Elements", data=ed)
    h5_file.close()
