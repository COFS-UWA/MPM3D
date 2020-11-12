import h5py as h5

msh_filename = "rect_mesh"

msh_file = h5.File("../Asset/" + msh_filename + ".h5", "r")

msh_grp = msh_file["Mesh"]

elem_grp = msh_grp["Elements"]
node_grp = msh_grp["Nodes"]

node_num = len(node_grp)
node_elem_pairs = []
for n_id in range(node_num):
    node_elem_pairs.append([])

elem_num = len(elem_grp)
for e_id in range(elem_num):
    elem = elem_grp[e_id]
    # node_elem_pairs[elem["n1"]-1].append(3 * e_id)
    # node_elem_pairs[elem["n2"]-1].append(3 * e_id + 1)
    # node_elem_pairs[elem["n3"]-1].append(3 * e_id + 2)
    node_elem_pairs[elem["n1"]-1].append(e_id)
    node_elem_pairs[elem["n2"]-1].append(e_id)
    node_elem_pairs[elem["n3"]-1].append(e_id)

msh_file.close()

# output to text file
res_file = open("../Build/TestsParallel/node_elem_pairs.txt", "w")
for n_id in range(node_num):
    res_file.write("%d: " % n_id)
    for e_id in node_elem_pairs[n_id]:
        res_file.write("%d, " % e_id)
    res_file.write("\n")
res_file.close()
