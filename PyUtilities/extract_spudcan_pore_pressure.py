import math
import sys
import h5py as py
import matplotlib.pyplot as plt
from multiprocessing import Pool

import find_closest_pcls as fp

file_name = "t3d_chm_mt_spudcan_-7"
pore_pos_offset = (0.0, 0.0, -0.3)

class AvgVarData:
    def __init__(self, pcl_ids, pcl_vols, pcl_xs, pcl_ys, pcl_zs, pcl_ps, \
                 pcl_num, closest_pcls_set, proc_num, proc_id):
        self.pcl_ids = pcl_ids
        self.pcl_vols = pcl_vols
        self.pcl_xs = pcl_xs
        self.pcl_ys = pcl_ys
        self.pcl_zs = pcl_zs
        self.pcl_ps = pcl_ps
        self.pcl_num = pcl_num
        self.closest_pcls_set = closest_pcls_set
        self.proc_num = proc_num
        self.proc_id = proc_id

def cal_avg_var(data):
    pcl_ids = data.pcl_ids
    pcl_vols = data.pcl_vols
    pcl_xs = data.pcl_xs
    pcl_ys = data.pcl_ys
    pcl_zs = data.pcl_zs
    pcl_ps = data.pcl_ps
    closest_pcls_set = data.closest_pcls_set
    all_vol = 0.0
    avg_x = 0.0
    avg_y = 0.0
    avg_z = 0.0
    avg_pore = 0.0
    #print("proc %d start %d." % (data.proc_id, data.pcl_num * data.proc_id // data.proc_num))
    #print("proc %d end %d." % (data.proc_id, data.pcl_num * (data.proc_id + 1) // data.proc_num))
    for i in range(data.pcl_num * data.proc_id // data.proc_num, \
                   data.pcl_num * (data.proc_id + 1) // data.proc_num):
        if pcl_ids[i] in closest_pcls_set:
            p_vol = pcl_vols[i]
            all_vol += p_vol
            avg_x += pcl_xs[i] * p_vol
            avg_y += pcl_ys[i] * p_vol
            avg_z += pcl_zs[i] * p_vol
            avg_pore += pcl_ps[i] * p_vol
    return ( all_vol, avg_x, avg_y, avg_z, avg_pore )

if __name__ == "__main__":
    # extract offset number spudcan tip
    tip_initial_pos = (0.0, 0.0, 0.0)
    closest_pcl_num = 10
    proc_num = 5

    data_file = open("../Build/TestsParallel/" + file_name + "_pore.csv", "w")
    data_file.write("pos_x, pos_y, pos_z, avg_x, avg_y, avg_z, avg_pore\n")

    # Numerical result
    hdf5_file = py.File("../Build/TestsParallel/" + file_name + ".h5", "r")

    pore_pos = (tip_initial_pos[0] + pore_pos_offset[0], \
                tip_initial_pos[1] + pore_pos_offset[1], \
                tip_initial_pos[2] + pore_pos_offset[2])
    
    # get closest points
    # md_p_dset = hdf5_file['ModelData']['ParticleData']['field']
    # closest_pcls, pcls_dist = fp.get_closest_pcls(md_p_dset, pore_pos, closest_pcl_num)
    # if len(closest_pcls) == 0:
        # raise UserWarning("Cannot get pcl in frame %d!\n" % th_id)
    # print(closest_pcls)

    closest_pcls = [651110, 651111, 651197, 651023, 651109, 651198, 651024, 651196, 651022, 651112]
    closest_pcls_set = set(closest_pcls)

    th_grp = hdf5_file['TimeHistory']['penetration']
    th_num = th_grp.attrs['output_num']

    rb_z = [0.0]
    pcl_pore = []
    ini_x = 0.0
    ini_y = 0.0
    ini_z = 0.0
    is_init = False
    proc_pool = Pool(processes = proc_num)
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
        proc_data = []
        for i in range(proc_num):
            proc_data.append(AvgVarData(pcl_ids, pcl_vols, pcl_xs, \
                                        pcl_ys, pcl_zs, pcl_ps, \
                                        pcl_num, closest_pcls_set, \
                                        proc_num, i))
        res = proc_pool.map(cal_avg_var, proc_data)
        # combine res
        all_vol = 0.0
        avg_x = 0.0
        avg_y = 0.0
        avg_z = 0.0
        avg_pore = 0.0
        # for i in range(pcl_num):
            # if pcl_ids[i] in closest_pcls:
                # p_vol = pcl_vols[i]
                # all_vol += p_vol
                # avg_x += pcl_xs[i] * p_vol
                # avg_y += pcl_ys[i] * p_vol
                # avg_z += pcl_zs[i] * p_vol
                # avg_pore += pcl_ps[i] * p_vol
            # if i % 100000 == (100000-1):
                # print("%d pcls processed." % (i + 1))
        for i in range(proc_num):
            proc_res = res[i];
            all_vol += proc_res[0]
            avg_x += proc_res[1]
            avg_y += proc_res[2]
            avg_z += proc_res[3]
            avg_pore += proc_res[4]
        if all_vol != 0.0:
            avg_x /= all_vol
            avg_y /= all_vol
            avg_z /= all_vol
            avg_pore /= all_vol
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
