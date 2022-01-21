import math
import h5py as py
import matplotlib.pyplot as plt

#file_name = "t3d_chm_mt_cap_compression"
file_name = "t3d_chm_tbb_cap_compression"

# Numerical result
hdf5_file = py.File("../Build/TestsParallel/" + file_name + ".h5", "r")

th_grp = hdf5_file['TimeHistory']['compression']
th_num = th_grp.attrs['output_num']

rb_z = [0.0]
rb_fz = [0.0]
is_init = False
ini_z = 0.0
for th_id in range(th_num):
    rb_grp = th_grp['frame_%d' % th_id]['RigidCylinder']
    cen_z = rb_grp.attrs['z']
    rf_z = rb_grp.attrs['fz_cont']
    if not is_init:
        ini_z = cen_z
        is_init = True
    else:
        rb_z.append(ini_z - cen_z)
        rb_fz.append(rf_z)

hdf5_file.close()

data_file = open("../Build/TestsParallel/" + file_name + "_rf.csv", "w")
for i in range(len(rb_z)):
    data_file.write("%f, %f\n" % (rb_z[i], rb_fz[i]))
data_file.close()

print(rb_z)
print(rb_fz)

fig = plt.figure()
plot1 = fig.subplots(1, 1)
line1, = plot1.plot(rb_z, rb_fz)

plt.xlim(rb_z[0], rb_z[-1])
plt.ylim(0.0)

plt.show()
