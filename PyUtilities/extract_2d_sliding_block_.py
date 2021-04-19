import h5py as py
import matplotlib.pyplot as plt

# Numerical result
hdf5_file = py.File("../Build/TestsParallel/t2d_me_mt_sliding_down_slope.h5", "r")

th_grp = hdf5_file['TimeHistory']['slide']
th_num = th_grp.attrs['output_num']

cal_time = []
rb_cfx = []
for th_id in range(th_num):
    frame_grp = th_grp['frame_%d' % th_id]
    cal_time.append(frame_grp.attrs['total_time'])
    rb_grp = frame_grp['RigidRect']
    rb_cfx.append(rb_grp.attrs['fx_contact'])

hdf5_file.close()

data_file = open("../Build/TestsParallel/t2d_me_mt_sliding_down_slope.csv", "w")
data_file.write("time, fx_contact,\n")
for i in range(len(cal_time)):
    data_file.write("%f, %f\n" % (cal_time[i], rb_cfx[i]))
data_file.close()

fig = plt.figure()
plot1 = fig.subplots(1, 1)
line1, = plot1.plot(cal_time, rb_cfx)

cal_time_range = [min(cal_time), max(cal_time)]
#rb_x_range = [min(rb_x), max(rb_x)]
plt.xlim(cal_time_range)
#plt.ylim(rb_x_range)

#plt.legend(handles=[line1,line2], labels=['MPM', 'Analytical'])
plt.show()
