import h5py as py
import matplotlib.pyplot as plt

# Numerical result
hdf5_file = py.File("../Build/TestsParallel/t3d_me_mt_cap_compression.h5", "r")

th_grp = hdf5_file['TimeHistory']['compression']
th_num = th_grp.attrs['output_num']

rb_z = []
rb_fz = []
is_init = False
ini_z = 0.0
for th_id in range(th_num):
    rb_grp = th_grp['frame_%d' % th_id]['RigidCylinder']
    cen_z = rb_grp.attrs['z']
    rf_z = rb_grp.attrs['fz']
    if not is_init:
        ini_z = cen_z
        is_init = True
    else:
        rb_z.append(ini_z - cen_z)
        rb_fz.append(rf_z)

hdf5_file.close()

data_file = open("../Build/TestsParallel/t3d_me_mt_cap_compression_rf.csv", "w")
for i in range(len(rb_z)):
    data_file.write("%f, %f\n" % (rb_z[i], rb_fz[i]))
data_file.close()

print(rb_z)
print(rb_fz)

fig = plt.figure()
plot1 = fig.subplots(1, 1)
line1, = plot1.plot(rb_fz, rb_z)

# Analytical solution
k_com = 1000.0 * 0.2 * 0.2 / 1.0
z_range = [min(rb_z), max(rb_z)]
fz_range = [z_range[0] * k_com, z_range[1] * k_com]
line2, = plot1.plot(fz_range, z_range)

plt.xlim(fz_range)
plt.ylim(z_range)

plt.legend(handles=[line1, line2], labels=['MPM', 'Analytical Solution'])
plt.show()
