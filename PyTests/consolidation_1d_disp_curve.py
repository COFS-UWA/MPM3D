import numpy as np
import matplotlib.pyplot as plt
import h5py as py

from OneDConsolidation import OneDConsolidation

fig = plt.figure()
plot1 = fig.subplots(1, 1)
plot1.set_xlabel("time")
plot1.set_ylabel("displacement")

out_time = []
pcl_var = []

# numerical solution
hdf5_file = py.File("..\\Build\\Tests\\t3d_chm_s_1d_consolidation.h5", "r")
th_grp = hdf5_file['TimeHistory']['consolidation']

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
    pcl_fld = pcl_dset[728]
    var = pcl_fld['z']
    if not is_init:
        init_z = var
        is_init = True
    var = var - init_z
    pcl_var.append(var)

hdf5_file.close()

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
time = 10.0 # time of consolidation
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

with open("consolidation_disp_ana.csv", "w") as out_file:
    for i in range(len(t_list)):
        out_file.write("%f, %f\n" % (t_list[i], u_list[i]))

plot1.plot(t_list, u_list, 'r--')
plt.legend(handles=[line1, line2], labels=['MPM', 'Analytical Solution'])
plt.show()
