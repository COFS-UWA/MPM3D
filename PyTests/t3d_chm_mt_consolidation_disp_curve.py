import numpy as np
import matplotlib.pyplot as plt
import h5py as py

from OneDConsolidation import OneDConsolidation

file_name = "t3d_chm_tbb_1d_consolidation_80kPa"

fig = plt.figure()
plot1 = fig.subplots(1, 1)
plot1.set_xlabel("time")
plot1.set_ylabel("displacement")

out_time = []
pcl_var = []

# numerical solution
hdf5_file = py.File("../Build/TestsParallel/" + file_name + ".h5", "r")
th_grp = hdf5_file['TimeHistory']['consolidation']

csv_file = open("../Build/TestsParallel/" + file_name + "_disp.csv", "w")

output_num = th_grp.attrs['output_num']
is_init = False
init_z = 0.0
for t_id in range(output_num):
    # frame
    frame_grp = th_grp['frame_%d' % t_id]
    frame_time = frame_grp.attrs['total_time']
    out_time.append(frame_time)
    # particle
    pcl_dset = frame_grp['ParticleData']['field']
    pcl_ids = pcl_dset['id']
    pcl_zs = pcl_dset['z']
    for p_id in range(len(pcl_dset)):
        if pcl_ids[p_id] == 2559:
            var = pcl_zs[p_id]
            if not is_init:
                init_z = var
                is_init = True
            var = var - init_z
            pcl_var.append(var)

hdf5_file.close()

csv_file.write("mpm,\n")
for o_id in range(output_num):
    csv_file.write("%f, %f\n" % (out_time[o_id], pcl_var[o_id]))

line1, = plot1.plot(out_time, pcl_var)

# analytical solution
u0 = 80.0e3
H = 1.0
E = 21.0e6
niu = 0.0 # possion ratio
kv = 5.0e-13
miu = 1.0e-3 # dynamic viscosity
time = 40.0 # time of consolidation

Es = (1 - niu) / (1 + niu) / (1 - 2.0*niu) * E # Es = (1-v) / (1 + v) / (1-2v) * E
Cv = kv * Es / miu
con_res = OneDConsolidation(Cv, Es, u0, H)
data_num = 100
t_list = np.zeros(data_num + 2)
u_list = np.zeros(data_num + 2)
t_list[0] = 0.0
u_list[0] = 0.0
t_list[1] = 0.0 # time for equilibrium
u_list[1] = u_list[0]
for i in range(data_num):
    t_list[i + 2] = time * float(i) / float(data_num)
    u_list[i + 2] = con_res.calSettlement(t_list[i + 2])
    t_list[i + 2] += t_list[1]

print("Cv: %f" % Cv)

csv_file.write("analytical,\n")
for o_id in range(data_num):
    csv_file.write("%f, %f\n" % (t_list[o_id], u_list[o_id]))

line2, = plot1.plot(t_list, u_list, 'r--')
with open("consolidation_disp_ana_SE.csv", "w") as out_file:
    for i in range(len(t_list)):
        out_file.write("%f, %f\n" % (t_list[i], u_list[i]))

plt.legend(handles=[line1, line2], labels=['MPM', 'Analytical Solution'])
plt.show()
