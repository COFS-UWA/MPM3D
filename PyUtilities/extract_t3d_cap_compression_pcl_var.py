import h5py as py
import matplotlib.pyplot as plt

# Numerical result
hdf5_file = py.File("../Build/TestsParallel/t3d_me_mt_triaxial_compression.h5", "r")

th_grp = hdf5_file['TimeHistory']['compression']
th_num = th_grp.attrs['output_num']

pcl_ids = [1985, 2307, 2369]

rb_z = []
rb_fz = []
is_init = False
ini_z = 0.0
pcl_e11 = []
pcl_e33 = []
pcl_s11 = []
pcl_s33 = []
pcl_pi = []
for p_id in range(len(pcl_ids)):
    pcl_e11.append([])
    pcl_e33.append([])
    pcl_s11.append([])
    pcl_s33.append([])
    pcl_pi.append([])

for th_id in range(th_num):
    rb_grp = th_grp['frame_%d' % th_id]['RigidCylinder']
    cen_z = rb_grp.attrs['z']
    rf_z = rb_grp.attrs['fz']
    if not is_init:
        ini_z = cen_z
        is_init = True
    rb_z.append(ini_z - cen_z)
    rb_fz.append(rf_z)
    pcl_grp = th_grp['frame_%d' % th_id]['ParticleData']
    pcl_fld = pcl_grp['field']
    mat_grp = th_grp['frame_%d' % th_id]['MaterialModel']
    mat_fld = mat_grp['SandHypoplasticityStb']
    for p_id in range(len(pcl_ids)):
        pcl_var = pcl_fld[pcl_ids[p_id]]
        pcl_e11[p_id].append(pcl_var['e11'])
        pcl_e33[p_id].append(-pcl_var['e33'])
        pcl_s11[p_id].append(pcl_var['s11'])
        pcl_s33[p_id].append(pcl_var['s33'])
        mat_id = pcl_var['mat_id']
        mat_var = mat_fld[mat_id]
        #pcl_pi[p_id].append(mat_var['p_i'])

hdf5_file.close()

# data_file = open("../Build/TestsParallel/t3d_me_mt_cap_compression_e33.csv", "w")
# for i in range(len(rb_z)):
    # data_file.write("%f, %f, %f\n" % (rb_z[i], pcl_e33[i], pcl_s33[i]))
# data_file.close()

fig = plt.figure()
plot1 = fig.subplots(1, 1)

lines = []
line_lb = []

# Analytical solution
#line1, = plot1.plot(rb_z, rb_z)
#lines.append(line1)
#line_lb.append('Acc')

for p_id in range(len(pcl_ids)):
    line2, = plot1.plot(rb_z, pcl_e11[p_id])
    #line2, = plot1.plot(rb_z, pcl_e33[p_id])
    #line2, = plot1.plot(rb_z, pcl_s11[p_id])
    line2, = plot1.plot(rb_z, pcl_s33[p_id])
    #line2, = plot1.plot(rb_z, pcl_pi[p_id])
    lines.append(line2)
    line_lb.append('%d' % pcl_ids[p_id])

rb_z_range = [min(rb_z), max(rb_z)]
plt.xlim(rb_z_range)
#plt.ylim(rb_z_range)

plt.legend(handles=lines, labels=line_lb)
plt.show()
