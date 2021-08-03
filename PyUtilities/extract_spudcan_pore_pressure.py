import math
import sys
import h5py as py
import matplotlib.pyplot as plt
import find_closest_pcls as fp

# extract offset number spudcan tip
pore_pos_offset = (0.0, 0.0, -0.1)
tip_initial_pos = (0.0, 0.0, 0.0)
closest_pcl_num = 10

data_file = open("../Build/TestsParallel/t3d_chm_mt_spudcan_pore.csv", "w")
data_file.write("pos_x, pos_y, pos_z, avg_x, avg_y, avg_z, avg_pore\n")

# Numerical result
hdf5_file = py.File("../Build/TestsParallel/t3d_chm_mt_spudcan.h5", "r")

pore_pos = (tip_initial_pos[0] + pore_pos_offset[0], \
            tip_initial_pos[1] + pore_pos_offset[1], \
            tip_initial_pos[2] + pore_pos_offset[2])

md_p_dset = hdf5_file['ModelData']['ParticleData']['field']
closest_pcls, pcls_dist = fp.get_closest_pcls(md_p_dset, pore_pos, closest_pcl_num)
if len(closest_pcls) == 0:
    raise UserWarning("Cannot get pcl in frame %d!\n" % th_id)
print(closest_pcls)
closest_pcls = [651110, 651111, 651197, 651023, 651109, 651198, 651024, 651196, 651022, 651112]

th_grp = hdf5_file['TimeHistory']['penetration']
th_num = th_grp.attrs['output_num']

rb_z = [0.0]
pcl_pore = []
ini_x = 0.0
ini_y = 0.0
ini_z = 0.0
is_init = False
for th_id in range(th_num):
    frame_grp = th_grp['frame_%d' % th_id]
    rb_grp = frame_grp['RigidObjectByT3DMesh']
    cen_x = rb_grp.attrs['x']
    cen_y = rb_grp.attrs['y']
    cen_z = rb_grp.attrs['z']
    if not is_init:
        ini_x = cen_x
        ini_y = cen_y
        ini_z = cen_z
        is_init = True
    else:
        rb_z.append(ini_z - cen_z)
    #
    pcl_dset = frame_grp['ParticleData']['field']
    pcl_ids = pcl_dset[:]['id']
    pcl_xs = pcl_dset[:]['x']
    pcl_ys = pcl_dset[:]['y']
    pcl_zs = pcl_dset[:]['z']
    pcl_vols = pcl_dset[:]['vol']
    pcl_ps = pcl_dset[:]['p']
    pcl_num = len(pcl_dset)
    all_vol = 0.0
    avg_x = 0.0
    avg_y = 0.0
    avg_z = 0.0
    avg_pore = 0.0
    for i in range(pcl_num):
        if pcl_ids[i] in closest_pcls:
            p_vol = pcl_vols[i]
            all_vol += p_vol
            avg_x += pcl_xs[i] * p_vol
            avg_y += pcl_ys[i] * p_vol
            avg_z += pcl_zs[i] * p_vol
            avg_pore += pcl_ps[i] * p_vol
        if i % 100000 == (100000-1):
            print("%d pcls processed." % (i + 1))
    avg_pore /= all_vol
    avg_x /= all_vol
    avg_y /= all_vol
    avg_z /= all_vol
    pcl_pore.append(avg_pore)
    # output data
    data_file.write("%f, %f, %f, %f, %f, %f, %f\n" % \
        (pore_pos[0] + cen_x - ini_x, pore_pos[1] + cen_y - ini_y, \
         pore_pos[2] + cen_z - ini_z, avg_x, avg_y, avg_z, avg_pore))
    print("Complete frame %d." % th_id)

hdf5_file.close()
data_file.close()

# plot data
fig = plt.figure()
plot1 = fig.subplots(1, 1)
line1, = plot1.plot(pcl_pore, rb_z)
plt.ylim(rb_z[0], rb_z[-1])
plt.show()
