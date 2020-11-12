import h5py as h5

msh_filename = "brick_mesh_1.00_2x2x10"

msh_file = h5.File("../Asset/" + msh_filename + ".h5", "r")
res_file = open("../Build/TestsParallel/t3d_mt_model_py.txt", "w")

msh_grp = msh_file["Mesh"]

elem_grp = msh_grp["Elements"]
node_grp = msh_grp["Nodes"]

node_num = len(node_grp)

for n_id in range(node_num):
    node = node_grp[n_id]
    res_file.write("%d, %f, %f, %f,\n" % (n_id, node["x"], node["y"], node["z"]))

elem_num = len(elem_grp)

for e_id in range(elem_num):
    elem = elem_grp[e_id]
    res_file.write("%d, %d, %d, %d, %d,\n" % (e_id, elem["n1"], elem["n2"], elem["n3"], elem["n4"]))

elem_pairs = []
node_elem_pairs = []
for n_id in range(node_num):
    elem_pairs.append([])
    node_elem_pairs.append([])

elem_num = len(elem_grp)
for e_id in range(elem_num):
    elem = elem_grp[e_id]
    elem_pairs[elem["n1"]].append(e_id)
    elem_pairs[elem["n2"]].append(e_id)
    elem_pairs[elem["n3"]].append(e_id)
    elem_pairs[elem["n4"]].append(e_id)
    node_elem_pairs[elem["n1"]].append(4 * e_id)
    node_elem_pairs[elem["n2"]].append(4 * e_id + 1)
    node_elem_pairs[elem["n3"]].append(4 * e_id + 2)
    node_elem_pairs[elem["n4"]].append(4 * e_id + 3)

msh_file.close()

# output to text file
start_id = 0;
end_id = 0;
for n_id in range(node_num):
    end_id = start_id + len(elem_pairs[n_id]);
    res_file.write("%d, %d, %d:\n" % (n_id, start_id, end_id))
    for e_id in elem_pairs[n_id]:
        res_file.write("%d, " % e_id)
    res_file.write("\n")
    for ne_id in node_elem_pairs[n_id]:
        res_file.write("%d, " % ne_id)
    res_file.write("\n")
    start_id = end_id
res_file.close()
