import h5py as py
import matplotlib.pyplot as plt

file_name = "t3d_chm_mt_cap_compression"
#file_name = "t3d_chm_tbb_cap_compression"
#mm_name = 'SandHypoplasticityStb'
mm_name = 'Norsand'

# Numerical result
hdf5_file = py.File("../Build/TestsParallel/" + file_name + ".h5", "r")

th_grp = hdf5_file['TimeHistory']['compression']
th_num = th_grp.attrs['output_num']

is_init = False
ini_z = 0.0
rb_z = []
soil_e = []
for th_id in range(th_num):
    frame_grd = th_grp['frame_%d' % th_id]
    pcl_fld = frame_grd['ParticleData']['field']
    pcl_vols = pcl_fld['vol']
    pcl_mat_ids = pcl_fld['mat_id']
    mat_fld = frame_grd['MaterialModel'][mm_name]
    mat_es = mat_fld['e']
    pcl_tvv = 0.0
    pcl_tv = 0.0
    for p_id in range(len(pcl_fld)):
        pcl_vol = pcl_vols[p_id]
        mat_id = pcl_mat_ids[p_id]
        pcl_e = mat_es[mat_id]
        pcl_tvv += pcl_e / (1.0 + pcl_e) * pcl_vol
        pcl_tv += pcl_vol
    rb_grp = frame_grd['RigidCylinder']
    cen_z = rb_grp.attrs['z']
    if not is_init:
        is_init = True
        ini_z = cen_z
    rb_z.append(ini_z - cen_z)
    soil_e.append(pcl_tvv / (pcl_tv - pcl_tvv))

hdf5_file.close()

data_file = open("../Build/TestsParallel/" + file_name + "_e.csv", "w")
for i in range(len(rb_z)):
    data_file.write("%f, %f\n" % (rb_z[i], soil_e[i]))
data_file.close()

print(rb_z)
print(soil_e)

fig = plt.figure()
plot1 = fig.subplots(1, 1)
line1, = plot1.plot(rb_z, soil_e)

#plt.xlim(z_range)
#plt.ylim(fz_range)

plt.show()
