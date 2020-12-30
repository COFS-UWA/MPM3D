import numpy as np
import h5py as h5

# Reformat .msh file created by Gmsh into hdf5 file
# Only 2D triangluar mesh

class MeshParser:
    def __init__(self):
        self.node_num = 0
        self.node_indices = None # node index
        self.node_coords = None # node x, y, z
        self.elem_num = 0
        self.elems = None # id, n1_id, n2_id, n3_id
        # for parsing
        self.line = None
    
    def __def__(self):
        pass
    
    # main function
    def parse(self, file_name):
        msh_file = open(msh_file_name, "r")
        while True:
            self.line = msh_file.readline()
            if self.line == "":
                return
            elif self.line == "$Nodes\n":
                self.parse_nodes(msh_file)
            elif self.line == "$Elements\n":
                self.parse_elements(msh_file)
        msh_file.close()

    def parse_nodes(self, msh_file):
        # parse header
        self.line = msh_file.readline()
        line_data = list(map(lambda x: int(x), self.line.strip('\n').split(' ')))
        self.node_num = line_data[1]
        print("node num: %i" % self.node_num)
        self.node_indices = np.zeros(self.node_num, dtype = np.int64)
        self.node_coords = np.zeros([self.node_num, 3], dtype = np.float64)
        # parse node block
        block_num = line_data[0]
        node_id = 0
        for block_id in range(block_num):
            # block header
            self.line = msh_file.readline()
            line_data = list(map(lambda x: int(x), self.line.strip('\n').split(' ')))
            block_node_num = line_data[3]
            # get node index
            for i in range(block_node_num):
                self.line = msh_file.readline()
                self.node_indices[node_id+i] = int(self.line.strip('\n'))
            # get node coords
            for i in range(block_node_num):
                self.line = msh_file.readline()
                node_coord = list(map(lambda x: float(x), self.line.strip('\n').split(' ')))
                self.node_coords[node_id+i][0] = node_coord[0]
                self.node_coords[node_id+i][1] = node_coord[1]
                self.node_coords[node_id+i][2] = node_coord[2]
            # finish
            node_id += block_node_num          
        # strip ender "$EndNodes"
        msh_file.readline()
    
    def parse_elements(self, msh_file):
        # parse header
        self.line = msh_file.readline()
        line_data = list(map(lambda x: int(x), self.line.strip('\n').split(' ')))
        # parse element block
        block_num = line_data[0]
        self.elem_num = 0
        self.elems = []
        for block_id in range(block_num):
            # block header
            self.line = msh_file.readline()
            line_data = list(map(lambda x: int(x), self.line.strip('\n').split(' ')))
            block_elem_type = line_data[2]
            block_elem_num = line_data[3]
            if block_elem_type != 2:
                # skip block
                for i in range(block_elem_num):
                    msh_file.readline()
            else:
                self.elem_num += block_elem_num
                for i in range(block_elem_num):
                    self.line = msh_file.readline()
                    elem_data = list(map(lambda x: int(x), self.line.strip('\n').strip(' ').split(' ')))
                    self.elems.append(elem_data)
        print("elem num: %i" % self.elem_num)
        # strip ender "$EndElements"
        msh_file.readline()



if __name__ == "__main__":
    #msh_filename = "strip_footing_soil_mesh"
    #msh_filename = "smaller_rect_mesh_denser"
    #msh_filename = "strip_footing_soil_mesh_denser"
    msh_filename = "rect_pipe_conference_mesh2"
    # parse
    msh_file_name = "../Asset/" + msh_filename + ".msh"
    msh_parser = MeshParser()
    msh_parser.parse(msh_file_name)
    # output to hdf5 file
    h5_file = h5.File("../Asset/" + msh_filename + ".h5", "w")
    mesh_grp = h5_file.create_group("Mesh")
    mesh_grp.attrs['type'] = "TriangleMesh"
    # write nodes
    node_ids = msh_parser.node_indices
    node_coords = msh_parser.node_coords
    node_type = np.dtype([
        ("index", np.int64),
        ("x", np.float64),
        ("y", np.float64)
        ])
    nd = np.empty(len(node_ids), dtype = node_type)
    for i in range(len(node_ids)):
        nd[i]["index"] = node_ids[i]
        nd[i]["x"] = node_coords[i][0]
        nd[i]["y"] = node_coords[i][1]
    n_dset = mesh_grp.create_dataset("Nodes", data = nd)
    # write elements
    elems = msh_parser.elems
    elem_type = np.dtype([
        ("index", np.int64),
        ("n1", np.int64),
        ("n2", np.int64),
        ("n3", np.int64)
        ])
    ed = np.empty(len(elems), dtype = elem_type)
    for i in range(len(elems)):
        ed[i]["index"] = elems[i][0]
        ed[i]["n1"] = elems[i][1]
        ed[i]["n2"] = elems[i][2]
        ed[i]["n3"] = elems[i][3]
    e_dset = mesh_grp.create_dataset("Elements", data=ed)
    h5_file.close()
