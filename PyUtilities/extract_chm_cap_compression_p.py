import h5py as py
import matplotlib.pyplot as plt

file_name = "t3d_chm_mt_cap_compression"
#file_name = "t3d_chm_tbb_cap_compression"

# Numerical result
hdf5_file = py.File("../Build/TestsParallel/" + file_name + ".h5", "r")

th_grp = hdf5_file['TimeHistory']['compression']
th_num = th_grp.attrs['output_num']

is_init = False
ini_z = 0.0
rb_z = []
pcl_p = []
for th_id in range(th_num):
    frame_grd = th_grp['frame_%d' % th_id]
    pcl_fld = frame_grd['ParticleData']['field']
    pcl_vols = pcl_fld['vol']
    pcl_pores = pcl_fld['p']
    pcl_tpv = 0.0
    pcl_tv = 0.0
    for p_id in range(len(pcl_fld)):
        pcl_vol = pcl_vols[p_id]
        pcl_pore = pcl_pores[p_id]
        pcl_tpv += pcl_pore * pcl_vol
        pcl_tv += pcl_vol
    rb_grp = frame_grd['RigidCylinder']
    cen_z = rb_grp.attrs['z']
    if not is_init:
        is_init = True
        ini_z = cen_z
    rb_z.append(ini_z - cen_z)
    pcl_p.append(pcl_tpv / pcl_tv)

hdf5_file.close()

data_file = open("../Build/TestsParallel/t3d_chm_mt_cap_compression_p.csv", "w")
for i in range(len(rb_z)):
    data_file.write("%f, %f\n" % (rb_z[i], pcl_p[i]))
data_file.close()

print(rb_z)
print(pcl_p)

fig = plt.figure()
plot1 = fig.subplots(1, 1)
line1, = plot1.plot(rb_z, pcl_p)

#plt.xlim(z_range)
#plt.ylim(fz_range)

plt.show()
