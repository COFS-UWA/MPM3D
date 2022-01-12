import numpy as np
import matplotlib.pyplot as plt
import h5py as py

from OneDConsolidation import OneDConsolidation

p_ids = [2, 82, 183, 253, 344, 418, 512, 643, 789]
frame_id = 80

fig = plt.figure()
plot1 = fig.subplots(1, 1)
plot1.set_xlabel("pore pressure")
plot1.set_ylabel("depth")

p_zs = []
p_ps = []

# numerical solution
csv_file = open("../Build/TestsParallel/t3d_chm_mt_1d_consolidation_pd_curve.csv", "w")
hdf5_file = py.File("../Build/TestsParallel/t3d_chm_mt_1d_consolidation.h5", "r")
th_grp = hdf5_file['TimeHistory']['consolidation']

output_num = th_grp.attrs['output_num']
frame_grp = th_grp['frame_%d' % frame_id]
frame_time = frame_grp.attrs['total_time']
csv_file.write("%f,\n" % frame_time)

pcl_dset = frame_grp['ParticleData']['field']
pcl_ids = pcl_dset['id']
pcl_zs = pcl_dset['z']
pcl_ps = pcl_dset['p']
pcl_num = len(pcl_dset)

for p_id in p_ids:
    for i in range(pcl_num):
        if p_id == pcl_ids[i]:
            p_zs.append(1.0 - pcl_zs[i])
            p_ps.append(pcl_ps[i])

print(p_zs)

hdf5_file.close()

csv_file.write("mpm,\nid, z, pore\n")
for i in range(len(p_ids)):
    csv_file.write("%f, %f, %f\n" % (p_ids[i], p_zs[i], p_ps[i]))

line1, = plot1.plot(p_ps, p_zs)

# analytical solution
u0 = 1.0
H = 1.0
E = 1000.0
niu = 0.0 # possion ratio
kv = 1.0e-4
miu = 1.0 # dynamic viscosity
data_num = 100

Es = (1 - niu) / (1 + niu) / (1 - 2.0*niu) * E # Es = (1-v) / (1 + v) / (1-2v) * E
Cv = kv * Es / miu
con_res = OneDConsolidation(Cv, Es, u0, H)
data_num += 1
z_list = np.zeros(data_num)
u_list = np.zeros(data_num)
for i in range(data_num):
    z_list[i] = H * float(i) / float(data_num)
    u_list[i] = con_res.calPorePressure(frame_time, z_list[i])

line2, = plot1.plot(u_list, z_list, 'r--')

plt.legend(handles=[line1, line2], labels=['MPM', 'Analytical Solution'])
plt.show()

csv_file.write("analytical,\n")
for o_id in range(data_num):
    csv_file.write("%f, %f\n" % (z_list[o_id], u_list[o_id]))
