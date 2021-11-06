import math
import h5py as py
import matplotlib.pyplot as plt

mm_name = "SandHypoplasticity"
hdf5_name = "../Build/TestsParallel/t3d_me_mt_cap_compression.h5"
pcl_id = 100

# Numerical result
hdf5_file = py.File(hdf5_name, "r")

th_grp = hdf5_file['TimeHistory']['compression']
th_num = th_grp.attrs['output_num']

ini_z = 0.0
is_init = False
rb_z = []
pcl_de11 = []
pcl_de22 = []
pcl_de33 = []
pcl_de12 = []
pcl_de23 = []
pcl_de31 = []
pcl_de_vol = []
pcl_de_dev = []
prev_pcl_e11 = 0.0
prev_pcl_e22 = 0.0
prev_pcl_e33 = 0.0
prev_pcl_e12 = 0.0
prev_pcl_e23 = 0.0
prev_pcl_e31 = 0.0
for th_id in range(th_num):
    rb_grp = th_grp['frame_%d' % th_id]['RigidCylinder']
    cen_z = rb_grp.attrs['z']
    if not is_init:
        is_init = True
        ini_z = cen_z
    rb_z.append(ini_z - cen_z)
    pcl_fld = th_grp['frame_%d' % th_id]['ParticleData']['field']
    pcl_fld_id = pcl_fld['id']
    pcl_fld_e11 = pcl_fld['e11']
    pcl_fld_e22 = pcl_fld['e22']
    pcl_fld_e33 = pcl_fld['e33']
    pcl_fld_e12 = pcl_fld['e12']
    pcl_fld_e23 = pcl_fld['e23']
    pcl_fld_e31 = pcl_fld['e31']
    pcl_num = len(pcl_fld_id)
    # mat_grp = th_grp['frame_%d' % th_id]['MaterialModel']
    # mat_fld = mat_grp[mm_name]
    for p_id in range(pcl_num):
        if pcl_id == pcl_fld_id[p_id]:
            p_de11 = pcl_fld_e11[p_id] - prev_pcl_e11
            p_de22 = pcl_fld_e22[p_id] - prev_pcl_e22
            p_de33 = pcl_fld_e33[p_id] - prev_pcl_e33
            p_de12 = pcl_fld_e12[p_id] - prev_pcl_e12
            p_de23 = pcl_fld_e23[p_id] - prev_pcl_e23
            p_de31 = pcl_fld_e31[p_id] - prev_pcl_e31
            prev_pcl_e11 = pcl_fld_e11[p_id]
            prev_pcl_e22 = pcl_fld_e22[p_id]
            prev_pcl_e33 = pcl_fld_e33[p_id]
            prev_pcl_e12 = pcl_fld_e12[p_id]
            prev_pcl_e23 = pcl_fld_e23[p_id]
            prev_pcl_e31 = pcl_fld_e31[p_id]
            pcl_de11.append(p_de11)
            pcl_de22.append(p_de22)
            pcl_de33.append(p_de33)
            pcl_de12.append(p_de12)
            pcl_de23.append(p_de23)
            pcl_de31.append(p_de31)
            pcl_de_vol.append(p_de11 + p_de22 + p_de33)
            pcl_de_dev.append(math.sqrt(0.5 * (p_de11*p_de11 + p_de22*p_de22 + p_de33*p_de33) \
                + 3.0 * (p_de12 * p_de12 + p_de23 * p_de23 + p_de31 * p_de31)))
            # mat_id = pcl_var['mat_id']
            # mat_var = mat_fld[mat_id]

hdf5_file.close()

# data_file = open("../Build/TestsParallel/t3d_me_mt_cap_compression_e33.csv", "w")
# for i in range(len(rb_z)):
    # data_file.write("%f, %f, %f\n" % (rb_z[i], pcl_e33[i], pcl_s33[i]))
# data_file.close()

fig = plt.figure()
plot1 = fig.subplots(1, 1)

line2, = plot1.plot(rb_z, pcl_de_vol)
#line2, = plot1.plot(rb_z, pcl_de_dev)

rb_z_range = [min(rb_z), max(rb_z)]
plt.xlim(rb_z_range)
plt.show()
