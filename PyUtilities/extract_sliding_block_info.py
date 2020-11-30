import h5py as py
import matplotlib.pyplot as plt

# Numerical result
hdf5_file = py.File("../Build/TestsParallel/t3d_me_mt_block_sliding.h5", "r")

th_grp = hdf5_file['TimeHistory']['sliding']
th_num = th_grp.attrs['output_num']

cal_time = []
rb_x = []
rb_y = []
rb_z = []
rb_cfx = []
rb_cfy = []
rb_cfz = []
for th_id in range(th_num):
    frame_grp = th_grp['frame_%d' % th_id]
    cal_time.append(frame_grp.attrs['current_time'])
    rb_grp = frame_grp['RigidCube']
    rb_x.append(rb_grp.attrs['x'])
    rb_y.append(rb_grp.attrs['y'])
    rb_z.append(rb_grp.attrs['z'])
    rb_cfx.append(rb_grp.attrs['fx_cont'])
    rb_cfy.append(rb_grp.attrs['fy_cont'])
    rb_cfz.append(rb_grp.attrs['fz_cont'])

hdf5_file.close()

data_file = open("../Build/TestsParallel/t3d_me_mt_block_sliding_res.csv", "w")
data_file.write("time, x, y, z, fx_cont, fy_cont, fz_cont,\n")
for i in range(len(cal_time)):
    data_file.write("%f, %f, %f, %f, %f, %f, %f,\n" % (cal_time[i], \
                    rb_x[i], rb_y[i], rb_z[i], rb_cfx[i], rb_cfy[i], rb_cfz[i]))
data_file.close()

fig = plt.figure()
plot1 = fig.subplots(1, 1)
line1, = plot1.plot(cal_time, rb_x)
# line1, = plot1.plot(cal_time, rb_cfz)

# Analytical solution
# a = 0.33333333333 # snooth
a = (3.0 - 3.0 * 0.2) / 9.0
an_rb_x = []
for i in range(len(cal_time)):
    an_rb_x.append(0.5 * a * cal_time[i] * cal_time[i] + rb_x[0])
line2, = plot1.plot(cal_time, an_rb_x, '--')

cal_time_range = [min(cal_time), max(cal_time)]
rb_x_range = [min(rb_x), max(rb_x)]
plt.xlim(cal_time_range)
plt.ylim(rb_x_range)

plt.legend(handles=[line1,line2], labels=['MPM', 'Analytical'])
plt.show()
