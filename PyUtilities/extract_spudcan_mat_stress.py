import math
import h5py as py
import matplotlib.pyplot as plt

file_name = "t3d_chm_mt_spudcan_cy_geo"

mat_id = 2542
spudcan_diameter = 3.0

# Numerical result
hdf5_file = py.File("../Build/TestsParallel/" + file_name + ".h5", "r")

th_grp = hdf5_file['TimeHistory']['geostatic']
th_num = th_grp.attrs['output_num']

frm_grp = th_grp['frame_20']
mat_grp = frm_grp['MaterialModel']
nor_dset = mat_grp['Norsand']
mat_s = nor_dset['s33']
for i in range(len(mat_s)):
    if mat_s[i] < -1000000.0:
        print("    %d," % i)

stress = []
dstrain = []
pos = []
void = []
e11_prev = 0.0
e22_prev = 0.0
e33_prev = 0.0
e12_prev = 0.0
e23_prev = 0.0
e31_prev = 0.0
spudcan_area = math.pi * spudcan_diameter * spudcan_diameter * 0.25
for th_id in range(th_num):
    frm_grp = th_grp['frame_%d' % th_id]
    pcl_grp = frm_grp['ParticleData']
    fld_dset = pcl_grp['field']
    mat_ids = fld_dset['mat_id']
    for i in range(len(mat_ids)):
        if mat_ids[i] == mat_id:
            pcl_dat = fld_dset[i]
            pcl_e11 = pcl_dat['e11']
            pcl_e22 = pcl_dat['e22']
            pcl_e33 = pcl_dat['e33']
            pcl_e12 = pcl_dat['e12']
            pcl_e23 = pcl_dat['e23']
            pcl_e31 = pcl_dat['e31']
            dstrain.append((pcl_e11 - e11_prev, pcl_e22 - e22_prev, pcl_e33 - e33_prev,
                            pcl_e12 - e12_prev, pcl_e23 - e23_prev, pcl_e31 - e31_prev))
            pos.append((pcl_dat['x'], pcl_dat['y'], pcl_dat['z']))
            e11_prev = pcl_e11
            e22_prev = pcl_e22
            e33_prev = pcl_e33
            e12_prev = pcl_e12
            e23_prev = pcl_e23
            e31_prev = pcl_e31
            print(th_id, i)
            break
    mat_grp = frm_grp['MaterialModel']
    nor_dset = mat_grp['Norsand'][mat_id]
    stress.append((nor_dset['s11'], nor_dset['s22'], nor_dset['s33'],
                   nor_dset['s12'], nor_dset['s23'], nor_dset['s31']))
    void.append(nor_dset['e'])

hdf5_file.close()

print("stress:\n")
print(stress)
print(len(stress))
print(dstrain)
print(len(dstrain))

data_file = open("../Build/TestsParallel/" + file_name + "_mat_stress.csv", "w")
data_file.write("s11, s22, s33, s12, s23, s31, e11, e22, e33, e12, e23, e31, x, y, z, e,\n")
for i in range(len(dstrain)):
    s_tmp = stress[i]
    e_tmp = dstrain[i]
    p = pos[i]
    data_file.write("%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f,\n" %
        (s_tmp[0], s_tmp[1], s_tmp[2], s_tmp[3], s_tmp[4], s_tmp[5],
         e_tmp[0], e_tmp[1], e_tmp[2], e_tmp[3], e_tmp[4], e_tmp[5],
         p[0], p[1], p[2], void[i]))
data_file.close()
