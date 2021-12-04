import numpy as np
import h5py

def generate_rect_mesh_nbym(row_num, col_num, hx, hy):
    # generate nodes
    nodes = []
    nr_num = row_num + 1
    nc_num = col_num + 1
    n_id = 0
    for r_id in range(nr_num):
        for c_id in range(nc_num):
            nodes.append([n_id, float(c_id)*hx, float(r_id)*hy])
            n_id += 1
    # generate elems
    elems = []
    e_id = 0
    for r_id in range(row_num):
        for c_id in range(col_num):
            n1_id = c_id + r_id * nc_num
            n2_id = n1_id + 1
            n3_id = c_id + (r_id + 1) * nc_num
            n4_id = n3_id + 1
            elems.append([e_id*2, n1_id, n2_id, n4_id])
            elems.append([e_id*2+1, n1_id, n4_id, n3_id])
            e_id += 1
    return nodes, elems
    
if __name__ == "__main__":
    # 1
    row_num = 1
    col_num = 4
    hx = 1.0
    hy = 1.0
    # 2
    # row_num = 12
    # col_num = 10
    # hx = 0.5
    # hy = 0.5
    # 3
    # row_num = 24
    # col_num = 20
    # hx = 0.25
    # hy = 0.25
    nodes, elems = generate_rect_mesh_nbym(row_num, col_num, hx, hy)
    h5_file = h5py.File("../Asset/rect_mesh_%dby%d.h5" % (row_num, col_num), "w")
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
