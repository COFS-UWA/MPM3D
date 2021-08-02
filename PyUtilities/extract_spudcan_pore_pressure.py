import math
import sys
import h5py as py
import matplotlib.pyplot as plt

# extract offset number spudcan tip
pore_pos_offset = [0.0, 0.0, -0.1]
tip_initial_pos = [0.0, 0.0, 0.0]

closest_pcl_num = 10

data_file = open("../Build/TestsParallel/t3d_chm_mt_spudcan_pore.csv", "w")
data_file.write("pos_x, pos_y, pos_z, avg_x, avg_y, avg_z, avg_pore\n")

# Numerical result
hdf5_file = py.File("../Build/TestsParallel/t3d_chm_mt_spudcan.h5", "r")

th_grp = hdf5_file['TimeHistory']['penetration']
th_num = th_grp.attrs['output_num']

rb_z = [0.0]
pcl_pore = [0.0]
ini_x = 0.0
ini_y = 0.0
ini_z = 0.0
is_init = False
for th_id in range(th_num):
    frame_grp = th_grp['frame_%d' % th_id]
    rb_grp = frame_grp['RigidBodyByT3DMesh']
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
    pore_pos = [tip_initial_pos[0] + cen_x - ini_x + pore_pos_offset[0], \
                tip_initial_pos[1] + cen_y - ini_y + pore_pos_offset[1], \
                tip_initial_pos[2] + cen_z - ini_z + pore_pos_offset[2]]
    pcl_dset = frame_grp['ParticleData']
    # get pcl close to pore_pos
    closest_pcls = []
    closest_pcl_dist2 = []
    min_pcl_dist2 = sys.float_info.max
    for p_id in range(len(pcl_dset)):
        p_data = pcl_dset[i]
        p_x = p_data['x']
        p_y = p_data['y']
        p_z = p_data['z']
        p_dist2 = (p_x - pore_pos[0]) * (p_x - pore_pos[0]) \
                + (p_y - pore_pos[1]) * (p_y - pore_pos[1]) \
                + (p_z - pore_pos[2]) * (p_z - pore_pos[2])
        if min_pcl_dist2 > p_dist2:
            min_pcl_dist2 = p_dist2
            if len(cloest_pcl) == 0:
                cloest_pcl.append(p_id)
                closest_pcl_dist2.append(min_pcl_dist2)
            else:
                for cp_id in range(len(cloest_pcl)):
                    if p_dist2 < closest_pcl_dist2[cp_id]:
                        cloest_pcl[cp_id].insert(cp_id, p_id)
                        closest_pcl_dist2[cp_id].insert(cp_id, p_dist2)
                        break
                if len(cloest_pcl) > closest_pcl_num:
                    cloest_pcl[cp_id].pop()
                    closest_pcl_dist2[cp_id].pop()
    if len(cloest_pcl) == 0:
        raise UserWarning("Cannot get pcl in frame %d!\n" % th_id)
    #
    all_vol = 0.0
    avg_pore = 0.0
    avg_x = 0.0
    avg_y = 0.0
    avg_z = 0.0
    for p_id in cloest_pcl:
        p_data = pcl_dset[p_id];
        p_vol = p_data['vol']
        all_vol += p_vol
        avg_pore += p_data['p'] * p_vol
        avg_x += p_data['x'] * p_vol
        avg_y += p_data['y'] * p_vol
        avg_z += p_data['z'] * p_vol
    avg_pore /= all_vol
    avg_x /= all_vol
    avg_y /= all_vol
    avg_z /= all_vol
    pcl_pore.append(avg_pore)
    # output data
    data_file.write("%f, %f, %f, %f, %f, %f, %f\n" % \
        (pore_pos[0], pore_pos[0], pore_pos[0], avg_x[0], avg_y[0], avg_z[0], avg_pore))
    print("Complete processing frame %d.\n" % th_id)

hdf5_file.close()
data_file.close()

# plot data
fig = plt.figure()
plot1 = fig.subplots(1, 1)
line1, = plot1.plot(pcl_pore, rb_z)
plt.ylim(rb_z[0], rb_z[-1])
plt.show()
