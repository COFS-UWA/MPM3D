import math as mh
import numpy as np
import h5py

def has_keyword(text_line):
    return text_line[0] == '*' and text_line[1] != '*'

# skip lines until the line contain keyword
def skip_lines_until_keyword(inp_file):
    while True:
        text_line = inp_file.readline()
        if text_line == "":
            return ""
        # keyword line
        if has_keyword(text_line):
            return text_line

def get_keyword(text_line):
    line_data = text_line.split(',')
    keyword_name = line_data[0][1:].strip()
    keyword_props = {}
    # if this keyword has attributes
    if len(line_data) > 1:
        for i in range(1, len(line_data)):
            prop_data = line_data[i].split('=')
            keyword_props[prop_data[0].strip()] = prop_data[1].strip()
    return keyword_name, keyword_props

def go_to_keyword(key_name, inp_file, current_line = None):
    keyword_name = ""
    keyword_props = {}
    if current_line and has_keyword(current_line):
        text_line = current_line
    else:
        text_line = skip_lines_until_keyword(inp_file)
    while text_line != "":
        keyword_name, keyword_props = get_keyword(text_line)
        # if (key_name == "Node"):
            # print(keyword_name)
            # print(keyword_name == key_name)
        if keyword_name == key_name:
            break
        text_line = skip_lines_until_keyword(inp_file)
    return keyword_name, keyword_props

def parse_nodes(inp_file):
    nodes = []
    text_line = inp_file.readline()
    while text_line != "" and text_line[0] != "*":
        node_data = text_line.split(',')
        nodes.append([int(node_data[0])-1, float(node_data[1]), float(node_data[2]), float(node_data[3])])
        text_line = inp_file.readline()
    return nodes, text_line

def parse_elements(inp_file):
    elems = []
    text_line = inp_file.readline()
    while text_line != "" and text_line[0] != "*":
        elem_data = text_line.split(',')
        elems.append([int(elem_data[0])-1, int(elem_data[1])-1, int(elem_data[2])-1, int(elem_data[3])-1, int(elem_data[4])-1])
        text_line = inp_file.readline()
    return elems, text_line

def get_mesh_from_inp(file_name, part_name):
    inp_file = open(file_name, 'r')
    # go the designated part
    while True:
        key_name, key_props = go_to_keyword("Part", inp_file)
        if key_name == "":
            break
        if key_props['name'] == part_name:
            break
    kn, kps = go_to_keyword("Node", inp_file)
    nodes, text_line = parse_nodes(inp_file)
    kn, kps = go_to_keyword("Element", inp_file, text_line)
    if (kps["type"] != "C3D4"):
        raise Exception("not C3D4 element.")
    elems, text_line = parse_elements(inp_file)
    inp_file.close()
    return nodes, elems

def output_mesh_to_hdf5(nodes, elems, file_name):
    h5_file = h5py.File(file_name, "w")
    mesh_grp = h5_file.create_group("Mesh")
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

def move_mesh_in_z_direction(nodes, z_off):
    for node in nodes:
        node[3] += z_off

def move_mesh_in_y_direction(nodes, y_off):
    for node in nodes:
        node[2] += y_off

def rotate_axses_by_angle(x_ang, y_ang, z_ang):
    x_ang = mh.radians(x_ang)
    y_ang = mh.radians(y_ang)
    z_ang = mh.radians(z_ang)
	#
    Kx = 0.0
    Ky = 0.0
    Kz = 0.0
    theta = mh.sqrt(x_ang*x_ang + y_ang*y_ang + z_ang*z_ang)
    if (theta != 0.0):
        Kx = x_ang / theta
        Ky = y_ang / theta
        Kz = z_ang / theta
    #
    q0 = mh.cos(0.5 * theta)
    q1 = Kx * mh.sin(0.5 * theta)
    q2 = Ky * mh.sin(0.5 * theta)
    q3 = Kz * mh.sin(0.5 * theta)
    #
    ix = np.zeros(3)
    iy = np.zeros(3)
    iz = np.zeros(3)
    ix[0] = 1.0 - 2.0 * (q2 * q2 + q3 * q3)
    ix[1] = 2.0 * (q1 * q2 - q0 * q3)
    ix[2] = 2.0 * (q1 * q3 + q0 * q2)
    iy[0] = 2.0 * (q1 * q2 + q0 * q3)
    iy[1] = 1.0 - 2.0 * (q3 * q3 + q1 * q1)
    iy[2] = 2.0 * (q2 * q3 - q0 * q1)
    iz[0] = 2.0 * (q1 * q3 - q0 * q2)
    iz[1] = 2.0 * (q2 * q3 + q0 * q1)
    iz[2] = 1.0 - 2.0 * (q1 * q1 + q2 * q2)
    return ix, iy, iz

def rotate_mesh(nodes, x_ang, y_ang, z_ang):
    ix, iy, iz = rotate_axses_by_angle(x_ang, y_ang, z_ang)
    for node in nodes:
        node_x_tmp = node[1]
        node_y_tmp = node[2]
        node_z_tmp = node[3]
        node[1] = ix[0] * node_x_tmp + ix[1] * node_y_tmp + ix[2] * node_z_tmp
        node[2] = iy[0] * node_x_tmp + iy[1] * node_y_tmp + iy[2] * node_z_tmp
        node[3] = iz[0] * node_x_tmp + iz[1] * node_y_tmp + iz[2] * node_z_tmp

if __name__ == "__main__":
    # text_line = "*Part, name=Part-1"
    # text_line = "*Node"
    # text_line = "*Element, type=C3D8R"
    # key_name, key_props = get_keyword(text_line)
    # print(key_name)
    # print(key_props)
    
    # file = open("Job-1.inp","r")
    # keyword_name, keyword_props = go_to_keyword("Element", file)
    # print(keyword_name)
    # print(keyword_props)
    # file.close()
    
    #nodes, elems = get_mesh_from_inp("../Asset/spudcan_model_flat_tip.inp", "Part-1")
    #output_mesh_to_hdf5(nodes, elems, "../Asset/spudcan_model_flat_tip.h5")
    # nodes, elems = get_mesh_from_inp("../Asset/spudcan_soil_quarter2.inp", "Part-1")
    # output_mesh_to_hdf5(nodes, elems, "../Asset/spudcan_soil_quarter2.h5")
    #nodes, elems = get_mesh_from_inp("../Asset/piezofoundation_soil_quarter.inp", "Part-1")
    #output_mesh_to_hdf5(nodes, elems, "../Asset/piezofoundation_soil_quarter.h5")
    # nodes, elems = get_mesh_from_inp("../Asset/cylinder_2x2_quad_model.inp", "Part-1")
    # output_mesh_to_hdf5(nodes, elems, "../Asset/cylinder_2x2_quad_model.h5")
    # nodes, elems = get_mesh_from_inp("../Asset/spudcan_soil_quarter_cylinder_8D.inp", "Part-1")
    # move_mesh_in_z_direction(nodes, -12.0)
    # output_mesh_to_hdf5(nodes, elems, "../Asset/spudcan_soil_quarter_cylinder_8D.h5")
    nodes, elems = get_mesh_from_inp("../Asset/spudcan_soil_half_cylinder_8R.inp", "Part-1")
    move_mesh_in_z_direction(nodes, -12.0)
    output_mesh_to_hdf5(nodes, elems, "../Asset/spudcan_soil_half_cylinder_8R.h5")
    # nodes, elems = get_mesh_from_inp("../Asset/weird_cylinder.inp", "Part-1")
    # output_mesh_to_hdf5(nodes, elems, "../Asset/weird_cylinder.h5")
    # nodes, elems = get_mesh_from_inp("../Asset/weird_block.inp", "Part-1")
    # rotate_mesh(nodes, 90.0, 0.0, 0.0)
    # move_mesh_in_y_direction(nodes, 3.0)
    # output_mesh_to_hdf5(nodes, elems, "../Asset/weird_block.h5")
    # nodes, elems = get_mesh_from_inp("../Asset/spudcan_model_Hossain_2006.inp", "Part-1")
    # output_mesh_to_hdf5(nodes, elems, "../Asset/spudcan_model_Hossain_2006.h5")
    # nodes, elems = get_mesh_from_inp("../Asset/spudcan_soil_quarter_Hossain_4D.inp", "Part-1")
    # move_mesh_in_z_direction(nodes, -15.0)
    # output_mesh_to_hdf5(nodes, elems, "../Asset/spudcan_soil_quarter_Hossain_4D.h5")
    # nodes, elems = get_mesh_from_inp("../Asset/cylinder_cap.inp", "Part-1")
    # move_mesh_in_z_direction(nodes, 1.0)
    # output_mesh_to_hdf5(nodes, elems, "../Asset/cylinder_cap.h5")
