import numpy as np
import matplotlib.pyplot as plt
import h5py as py

from BarAxialVibration import BarAxialVibration

out_time = []
pcl_var = []

# numerical solution
hdf5_file = py.File("../Build/Tests/t2d_me_s_1d_compression.h5", "r")
th_grp = hdf5_file['TimeHistory']['compression']

output_num = th_grp.attrs['output_num']
is_init = False
init_y = 0.0
for t_id in range(output_num):
    # frame
    frame_grp = th_grp['frame_%d' % t_id]
    frame_time = frame_grp.attrs['total_time']
    out_time.append(frame_time)
    # particle
    pcl_dset = frame_grp['ParticleData']['field']
    #pcl_fld = pcl_dset[21] #21 mid
    pcl_fld = pcl_dset[117] #117 bottom
    var = -pcl_fld['s22']
    if not is_init:
        init_y = pcl_fld['y']
        print(init_y)
        is_init = True
    pcl_var.append(var)

hdf5_file.close()

# analytical solution
H = 1.0
p0 = 0.1
bf = 0.0
E = 1000.0
density = 10.0
t_len = 1.0 # time length
data_num = 200
# cal data
bav = BarAxialVibration(H, p0, bf, E, density)
data_num += 1
t_ana = np.zeros(data_num)
s22_ana = np.zeros(data_num)
t_inv = t_len / float(data_num)
for i in range(data_num):
    t_ana[i] = t_inv * i
    s22_ana[i] = bav.stress(init_y, t_ana[i])

# plot var - time curve
fig = plt.figure()

plot1 = fig.subplots(1, 1)
plot1.set_xlabel("time")
plot1.set_ylabel("stress")

line1, = plot1.plot(out_time, pcl_var)
line2, = plot1.plot(t_ana, s22_ana)

plt.legend(handles=[line1, line2], labels=['MPM', 'Analytical Solution'])
plt.show()

# output to csv file
csv_file = open('t2d_1d_compression_stress.csv', 'w')

csv_file.write("MPM,\n")
# time
csv_file.write("time, ")
for ot in out_time:
    csv_file.write("%f, " % ot)
csv_file.write("\n")
# pcl vars
csv_file.write("s22, ")
for pv in pcl_var:
    csv_file.write("%f, " % pv)
csv_file.write("\n")

csv_file.write("Analytical solution,\n")
# time
csv_file.write("time, ")
for ot in t_ana:
    csv_file.write("%f, " % ot)
csv_file.write("\n")
# pcl vars
csv_file.write("s22, ")
for pv in s22_ana:
    csv_file.write("%f, " % pv)
csv_file.write("\n")

csv_file.close()
