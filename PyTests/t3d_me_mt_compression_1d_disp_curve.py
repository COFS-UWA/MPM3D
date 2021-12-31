import numpy as np
import matplotlib.pyplot as plt
import h5py as py

from BarAxialVibration import BarAxialVibration

out_time = []
pcl_var = []

# numerical solution
hdf5_file = py.File("../Build/TestsParallel/t3d_me_mt_1d_compression.h5", "r")
th_grp = hdf5_file['TimeHistory']['compression']

csv_file = open("../Build/TestsParallel/t3d_me_mt_1d_compression_disp.csv", "w")

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
    for i in range(len(pcl_dset)):
        if (pcl_ids[i] == 2559):
            var = pcl_zs[i]
            if not is_init:
                init_z = var
                is_init = True
            var = init_z - var
            pcl_var.append(var)
            break

hdf5_file.close()

csv_file.write("mpm,\n")
for o_id in range(output_num):
    csv_file.write("%f, %f\n" % (out_time[o_id], pcl_var[o_id]))

# analytical solution
H = 1.0
p0 = 10.0
bf = 0.0
E = 1000.0
density = 10.0
t_len = 2.0 # time length
data_num = 200
# cal data
bav = BarAxialVibration(H, p0, bf, E, density)
data_num += 1
t_ana = np.zeros(data_num)
u_ana = np.zeros(data_num)
t_inv = t_len / float(data_num)
for i in range(data_num):
    t_ana[i] = t_inv * i
    u_ana[i] = bav.displacement(H, t_ana[i])

csv_file.write("analytical,\n")
for o_id in range(data_num):
    csv_file.write("%f, %f\n" % (t_ana[o_id], u_ana[o_id]))

# plot var - time curve
fig = plt.figure()

plot1 = fig.subplots(1, 1)
plot1.set_xlabel("time")
plot1.set_ylabel("displacement")

line1, = plot1.plot(out_time, pcl_var)
line2, = plot1.plot(t_ana, u_ana)

plt.legend(handles=[line1, line2], labels=['MPM', 'Analytical Solution'])
plt.show()
