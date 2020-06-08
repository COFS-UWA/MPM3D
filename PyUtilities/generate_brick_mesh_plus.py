import numpy as np
import h5py

def add_elems_in_cube(n0, n1, n2, n3, n4, n5, n6, n7, elems):
    e_id_off = len(elems)
    elems.append([e_id_off+0, n0, n1, n2, n5])
    elems.append([e_id_off+1, n0, n2, n3, n7])
    elems.append([e_id_off+2, n5, n4, n7, n0])
    elems.append([e_id_off+3, n5, n7, n6, n2])
    elems.append([e_id_off+4, n5, n7, n2, n0])

def generate_brick_mesh(nx, ny, nz, xl, xu, yl, yu, zl, zu):
    # generate nodes
    dx = (xu - xl) / float(nx)
    dy = (yu - yl) / float(ny)
    dz = (zu - zl) / float(nz)
    nx += 1
    ny += 1
    nz += 1
    nodes = []
    n_id = 0
    for z_id in range(nz):
        for y_id in range(ny):
            for x_id in range(nx):
                nodes.append([n_id, xl + float(x_id)*dx, yl + float(y_id)*dy, zl + float(z_id)*dz])
                n_id += 1
    # generate elems
    ex = nx - 1
    ey = ny - 1
    ez = nz - 1
    elems = []
    for z_id in range(ez):
        for y_id in range(ey):
            for x_id in range(ex):
                n0 = x_id + y_id * nx + z_id * nx*ny
                n1 = n0 + 1
                n2 = n1 + nx
                n3 = n0 + nx
                n4 = n0 + nx*ny
                n5 = n4 + 1
                n6 = n5 + nx
                n7 = n4 + nx
                add_elems_in_cube(n0, n1, n2, n3, n4, n5, n6, n7, elems)
    return nodes, elems
    
if __name__ == "__main__":
    h5_file = h5py.File("..\\Asset\\brick_mesh_plus.h5", "w")
    mesh_grp = h5_file.create_group("Mesh")
    nodes, elems = generate_brick_mesh(2, 2, 10, 0.0, 1.0/10*2, 0.0, 1.0/10*2, 0.0, 1.0)
    # write nodes
    node_type = np.dtype([("index", np.int64), ("x", np.float64), ("y", np.float64), ("z", np.float64)])
    nd = np.empty(len(nodes), dtype=node_type)
    for i in range(len(nodes)):
        nd[i]["index"] = nodes[i][0]
        nd[i]["x"] = nodes[i][1]
        nd[i]["y"] = nodes[i][2]
        nd[i]["z"] = nodes[i][3]
    n_dset = mesh_grp.create_dataset("Nodes", data=nd)
    # write elements
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
