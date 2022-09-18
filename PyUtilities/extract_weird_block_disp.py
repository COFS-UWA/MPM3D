import math
import h5py as py
import matplotlib.pyplot as plt

file_name = "t3d_me_mt_weird_block_sliding"

# Numerical result
hdf5_file = py.File("../Build/TestsParallel/" + file_name + ".h5", "r")

th_grp = hdf5_file['TimeHistory']['sliding']
th_num = th_grp.attrs['output_num']

cur_time = []
rb_x = []
rb_z = []
for th_id in range(th_num):
    frm_grp = th_grp['frame_%d' % th_id]
    cur_time.append(frm_grp.attrs['current_time'])
    rb_grp = frm_grp['RigidObjectByT3DMesh']
    rb_x.append(rb_grp.attrs['x'])
    rb_z.append(rb_grp.attrs['z'])

hdf5_file.close()

data_file = open("../Build/TestsParallel/" + file_name + "_disp.csv", "w")
for th_id in range(th_num):
    data_file.write("%f, %f, %f,\n" % (cur_time[th_id], rb_x[th_id], rb_z[th_id]))
data_file.close()

# fig = plt.figure()
# plot1 = fig.subplots(1, 1)
# line1, = plot1.plot(rb_fz_norm, rb_z_norm)
# plt.xlim(0.0)
# plt.ylim(rb_z_norm[0], rb_z_norm[-1])
# plt.show()
