import math
import h5py as py
import matplotlib.pyplot as plt

spudcan_diameter = 3.0

# Numerical result
hdf5_file = py.File("../Build/TestsParallel/t3d_chm_mt_spudcan.h5", "r")
#hdf5_file = py.File("../Build/TestsParallel/t3d_me_mt_spudcan.h5", "r")

th_grp = hdf5_file['TimeHistory']['penetration']
th_num = th_grp.attrs['output_num']

rb_z = [0.0]
rb_fz = [0.0]
rb_z_norm = [0.0]
rb_fz_norm = [0.0]
is_init = False
ini_z = 0.0
spudcan_area = math.pi * spudcan_diameter * spudcan_diameter * 0.25
for th_id in range(th_num):
    rb_grp = th_grp['frame_%d' % th_id]['RigidObjectByT3DMesh']
    cen_z = rb_grp.attrs['z']
    rf_z = rb_grp.attrs['fz_cont']
    if not is_init:
        ini_z = cen_z
        is_init = True
    else:
        rb_z.append(ini_z - cen_z)
        rb_fz.append(rf_z*4.0)
        rb_z_norm.append((ini_z - cen_z)/spudcan_diameter)
        rb_fz_norm.append(rf_z*4.0/spudcan_area/1000.0)

hdf5_file.close()

data_file = open("../Build/TestsParallel/t3d_chm_mt_spudcan_rf.csv", "w")
for i in range(len(rb_z)):
    data_file.write("%f, %f, %f, %f\n" % (rb_z[i], rb_fz[i], rb_z_norm[i], rb_fz_norm[i]))
data_file.close()

print(rb_z)
print(rb_fz)

fig = plt.figure()
plot1 = fig.subplots(1, 1)
line1, = plot1.plot(rb_fz_norm, rb_z_norm)
plt.xlim(0.0)
plt.ylim(rb_z_norm[0], rb_z_norm[-1])
plt.show()
