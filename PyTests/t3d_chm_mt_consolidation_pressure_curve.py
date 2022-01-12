import numpy as np
import matplotlib.pyplot as plt
import h5py as py

from OneDConsolidation import OneDConsolidation

fig = plt.figure()
plot1 = fig.subplots(1, 1)
plot1.set_xlabel("time")
plot1.set_ylabel("pore pressure")

out_time = []
pcl_var = []

# numerical solution
hdf5_file = py.File("../Build/TestsParallel/t3d_chm_mt_1d_consolidation.h5", "r")
th_grp = hdf5_file['TimeHistory']['consolidation']

csv_file = open("../Build/TestsParallel/t3d_chm_mt_1d_consolidation_p.csv", "w")

output_num = th_grp.attrs['output_num']
pcl_z = 0.0
is_init = False
for t_id in range(output_num):
    # frame
    frame_grp = th_grp['frame_%d' % t_id]
    frame_time = frame_grp.attrs['total_time']
    out_time.append(frame_time)
    # particle
    pcl_dset = frame_grp['ParticleData']['field']
    pcl_ids = pcl_dset['id']
    pcl_zs = pcl_dset['z']
    pcl_ps = pcl_dset['p']
    for p_id in range(len(pcl_dset)):
        if pcl_ids[p_id] == 1: # at bottom
            if not is_init:
                pcl_z = pcl_zs[p_id]
                is_init = True
            var = pcl_ps[p_id]
            pcl_var.append(var)

hdf5_file.close()

csv_file.write("mpm,\n")
for o_id in range(output_num):
    csv_file.write("%f, %f\n" % (out_time[o_id], pcl_var[o_id]))

line1, = plot1.plot(out_time, pcl_var)

# analytical solution
u0 = 1.0
H = 1.0
E = 1000.0
niu = 0.0 # possion ratio
kv = 1.0e-4
miu = 1.0 # dynamic viscosity

Es = (1 - niu) / (1 + niu) / (1 - 2.0*niu) * E # Es = (1-v) / (1 + v) / (1-2v) * E
Cv = kv * Es / miu
con_res = OneDConsolidation(Cv, Es, u0, H)
time = 15.0 # time of consolidation
data_num = 100
t_list = np.zeros(data_num + 2)
u_list = np.zeros(data_num + 2)
t_list[0] = 0.0
u_list[0] = 0.0
t_list[1] = 0.0 # time for equilibrium
u_list[1] = u_list[0]
for i in range(data_num):
    t_list[i + 2] = time * float(i) / float(data_num)
    u_list[i + 2] = con_res.calPorePressure(t_list[i + 2], 1.0 - pcl_z)
    t_list[i + 2] += t_list[1]

line2, = plot1.plot(t_list, u_list, 'r--')

plt.legend(handles=[line1, line2], labels=['MPM', 'Analytical Solution'])
plt.show()

csv_file.write("analytical,\n")
for o_id in range(data_num):
    csv_file.write("%f, %f\n" % (t_ana[o_id], u_ana[o_id]))
