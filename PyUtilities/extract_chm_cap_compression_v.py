import math
import h5py as py
import matplotlib.pyplot as plt

file_name = "t3d_chm_mt_cap_compression"
#file_name = "t3d_chm_tbb_cap_compression"

# Numerical result
hdf5_file = py.File("../Build/TestsParallel/" + file_name + ".h5", "r")

th_grp = hdf5_file['TimeHistory']['compression']
th_num = th_grp.attrs['output_num']

frm_t = [0.0]
rb_vz = [0.0]
for th_id in range(th_num):
    frm_grp = th_grp['frame_%d' % th_id]
    frm_t.append(frm_grp.attrs['current_time'])
    rb_grp = frm_grp['RigidCylinder']
    rb_vz.append(rb_grp.attrs['vz'])

hdf5_file.close()

data_file = open("../Build/TestsParallel/" + file_name + "_v.csv", "w")
for i in range(th_num):
    data_file.write("%f, %f\n" % (frm_t[i], rb_vz[i]))
data_file.close()

print(frm_t)
print(rb_vz)

fig = plt.figure()
plot1 = fig.subplots(1, 1)
line1, = plot1.plot(frm_t, rb_vz)

#plt.xlim(rb_z[0], rb_z[-1])
#plt.ylim(0.0)

plt.show()
