import h5py as py
import matplotlib.pyplot as plt

# Numerical result
hdf5_file = py.File("../Build/TestsParallel/t2d_me_mt_cap_compression.h5", "r")

th_grp = hdf5_file['TimeHistory']['loading']
th_num = th_grp.attrs['output_num']

rb_y = []
rb_fy = []
is_init = False
ini_y = 0.0
for th_id in range(th_num):
    rb_grp = None
    rf_y = None
    if 'RigidCylinder' in th_grp['frame_%d' % th_id].keys():
        rb_grp = th_grp['frame_%d' % th_id]['RigidRect']
        rf_y = rb_grp.attrs['fy_contact']
    else:
        rb_grp = th_grp['frame_%d' % th_id]['RigidObjectByT2DMesh']
        rf_y = rb_grp.attrs['fy_cont']
    #rb_grp = th_grp['frame_%d' % th_id]['RigidObjectByT2DMesh'] # ['RigidRect'] ['RigidObjectByT2DMesh']
    #rf_y = rb_grp.attrs['fy_cont'] # ['fy_contact'] ['fy_cont']
    cen_y = rb_grp.attrs['y']
    if not is_init:
        ini_y = cen_y
        is_init = True
    else:
        rb_y.append(ini_y - cen_y)
        rb_fy.append(rf_y)

hdf5_file.close()

data_file = open("../Build/TestsParallel/t2d_me_mt_cap_compression_rf.csv", "w")
for i in range(len(rb_y)):
    data_file.write("%f, %f\n" % (rb_y[i], rb_fy[i]))
data_file.close()

print(rb_y)
print(rb_fy)

fig = plt.figure()
plot1 = fig.subplots(1, 1)
line1, = plot1.plot(rb_fy, rb_y)

# Analytical solution
k_com = 1000.0 * 0.2 / 1.0
y_range = [min(rb_y), max(rb_y)]
fy_range = [y_range[0] * k_com, y_range[1] * k_com]
line2, = plot1.plot(fy_range, y_range)

plt.xlim(fy_range)
plt.ylim(y_range)

plt.legend(handles=[line1,line2], labels=['MPM', 'Analytical Solution'])
plt.show()
