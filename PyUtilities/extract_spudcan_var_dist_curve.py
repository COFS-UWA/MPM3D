import math
import sys
import h5py as py
import matplotlib.pyplot as plt
from multiprocessing import Process, Array

path_name = "E:/"
#path_name = "../Build/TestsParallel/"
sim_name = "t3d_chm_mt_spudcan_cy_ucav-400_V16000"
var_name = "s33"
th_name = "penetration"
rb_name = "RigidObjectByT3DMesh"
th_id = 101
proc_num = 8

# (r0, z0, dr, dz)
# absolute rb_position
# var_depth = -3.0
# var_pos = [
    # (0.1, var_depth, 0.2, 0.2),
    # (0.3, var_depth, 0.2, 0.2),
    # (0.5, var_depth, 0.2, 0.2),
    # (0.75, var_depth, 0.2, 0.2),
    # (1.0, var_depth, 0.2, 0.2),
    # (1.25, var_depth, 0.2, 0.2),
    # (1.5, var_depth, 0.2, 0.2),
    # (1.7, var_depth, 0.2, 0.2),
    # (2.0, var_depth, 0.2, 0.2),
    # (2.5, var_depth, 0.2, 0.2),
    # (3.0, var_depth, 0.2, 0.2),
    # (5.0, var_depth, 0.2, 0.2),
    # (7.0, var_depth, 0.2, 0.2),
    # (8.0, var_depth, 0.2, 0.2),
# ]

# vertical
var_pos = [
    (0.1, -0.1, 0.2, 0.2),
    (0.1, -0.3, 0.2, 0.2),
    (0.1, -0.5, 0.2, 0.2),
    (0.1, -0.6, 0.2, 0.2),
    (0.1, -0.65, 0.2, 0.2),
    (0.1, -0.75, 0.2, 0.2),
    (0.1, -1.0, 0.2, 0.2),
    (0.1, -1.25, 0.2, 0.2),
    (0.1, -1.5, 0.2, 0.2),
    (0.1, -1.7, 0.2, 0.2),
    (0.1, -2.0, 0.2, 0.2),
    (0.1, -2.5, 0.2, 0.2),
    (0.1, -3.0, 0.2, 0.2),
    (0.1, -4.0, 0.2, 0.2),
    (0.1, -5.0, 0.2, 0.2),
    (0.1, -6.0, 0.2, 0.2),
    (0.1, -7.0, 0.2, 0.2),
    (0.1, -8.0, 0.2, 0.2),
    (0.1, -9.0, 0.2, 0.2),
    (0.1, -10.0, 0.2, 0.2),
    (0.1, -11.0, 0.2, 0.2),
    (0.1, -11.45, 0.2, 0.2),
]

# around spudcan
# var_pos = [
    # (0.1, -0.1, 0.2, 0.2),
    # (0.157, -0.1, 0.2, 0.2),
    # (0.3, -0.062, 0.2, 0.2),
    # (0.5, -0.008, 0.2, 0.2),
    # (0.7, 0.046, 0.2, 0.2),
    # (0.9, 0.1, 0.2, 0.2),
    # (1.1, 0.153, 0.2, 0.2),
    # (1.3, 0.206, 0.2, 0.2),
    # (1.5, 0.26, 0.2, 0.2),
    # (1.7, 0.26, 0.2, 0.2),
    # (1.9, 0.26, 0.2, 0.2),
# ]

###########################################################################################
###########################################################################################
def get_avg_var(pcl_xs, pcl_ys, pcl_zs, pcl_vols, pcl_vars, pcl_num, \
                model_var_r, model_var_z, model_var_dr, model_var_dz, \
                model_cen_x, model_cen_y):
    avg_var = 0.0
    eff_pcl_vol = 0.0
    eff_pcl_num = 0
    avg_rmin2 = (model_var_r - model_var_dr) * (model_var_r - model_var_dr)
    avg_rmax2 = (model_var_r + model_var_dr) * (model_var_r + model_var_dr)
    for pcl_id in range(pcl_num):
        md_p_dx = pcl_xs[pcl_id] - model_cen_x
        md_p_dy = pcl_ys[pcl_id] - model_cen_y
        md_p_dr2 = md_p_dx * md_p_dx + md_p_dy * md_p_dy
        md_p_dz = abs(pcl_zs[pcl_id] - model_var_z)
        if md_p_dr2 <= avg_rmax2 and md_p_dr2 >= avg_rmin2 and md_p_dz <= model_var_dz:
            avg_var += pcl_vars[pcl_id] * pcl_vols[pcl_id]
            eff_pcl_vol += pcl_vols[pcl_id]
            eff_pcl_num += 1
    if (eff_pcl_vol != 0.0):
        avg_var /= eff_pcl_vol
    # average pore pressure, effective pcl vol, effective pcl num
    return (avg_var, eff_pcl_vol, eff_pcl_num)

class AvgVarParam:
    def __init__(self, path_name, sim_name, th_name, th_id, \
                 rb_name, var_name, rb_ini_z, var_pos, \
                 model_cen_x = 0.0, model_cen_y = 0.0):
        self.path_name = path_name
        self.sim_name = sim_name
        self.th_name = th_name
        self.th_id = th_id
        self.rb_name = rb_name
        self.var_name = var_name
        self.rb_ini_z = rb_ini_z
        self.var_pos = var_pos
        self.model_cen_x = model_cen_x
        self.model_cen_y = model_cen_y

def proc_get_th_avg_var(pos_id_range, param, pcl_r, pcl_z, pcl_var):
    hdf5_file = py.File(param.path_name + param.sim_name + ".h5", "r")
    th_grp = hdf5_file['TimeHistory'][param.th_name]
    frame_grp = th_grp['frame_%d' % param.th_id]
    rb_grp = frame_grp[param.rb_name]
    rb_dz = rb_grp.attrs['z'] - param.rb_ini_z
    # cal avg pore pressure for each time histories
    for pos_id in range(pos_id_range[0], pos_id_range[1]):
        #
        avg_r = param.var_pos[pos_id][0]
        avg_z = param.var_pos[pos_id][1] + rb_dz
        pcl_r[pos_id] = avg_r
        pcl_z[pos_id] = avg_z
        pcl_dset = frame_grp['ParticleData']['field']
        pcl_xs = pcl_dset[:]['x']
        pcl_ys = pcl_dset[:]['y']
        pcl_zs = pcl_dset[:]['z']
        pcl_vols = pcl_dset[:]['vol']
        pcl_vars = pcl_dset[:][param.var_name]
        pcl_num = len(pcl_dset)
        avg_var, eff_pcl_vol, eff_pcl_num = get_avg_var( \
            pcl_xs, pcl_ys, pcl_zs, pcl_vols, pcl_vars, pcl_num, \
            avg_r, avg_z, param.var_pos[pos_id][2], param.var_pos[pos_id][3], \
            param.model_cen_x, param.model_cen_y)
        pcl_var[pos_id] = avg_var
        # print msg
        if pos_id_range[0] == 0:
            print("processed %d%%, %d pcls nearby." % (int((pos_id+1)/(pos_id_range[1]-pos_id_range[0])*100.0), eff_pcl_num))
    hdf5_file.close()

if __name__ == "__main__":
    hdf5_file = py.File(path_name + sim_name + ".h5", "r")
    th_grp = hdf5_file['TimeHistory'][th_name]
    frame_grp = th_grp['frame_0']
    rb_grp = frame_grp[rb_name]
    rb_ini_z = rb_grp.attrs['z']
    hdf5_file.close()
    
    var_num = len(var_pos)
    proc_pcl_r = Array('d', [0.0] * var_num)
    proc_pcl_z = Array('d', [0.0] * var_num)
    proc_pcl_var = Array('d', [0.0] * var_num)
    procs = []
    proc_param = AvgVarParam(path_name, sim_name, th_name, th_id, rb_name, var_name, rb_ini_z, var_pos)
    for proc_id in range(proc_num):
        p = Process(target = proc_get_th_avg_var, \
                    args = ([int(proc_id*var_num/proc_num), int((proc_id+1)*var_num/proc_num)], \
                             proc_param, proc_pcl_r, proc_pcl_z, proc_pcl_var))
        procs.append(p)
        p.start()
    
    print("All procs started.")
    
    for proc_id in range(proc_num):
        procs[proc_id].join()
        print("proc %d completed." % proc_id)
    
    print("All procs completed.")

    pcl_r = []
    pcl_z = []
    pcl_var = []
    csv_file = open(path_name + sim_name + "_" + var_name + ".csv", "w")
    csv_file.write("pcl_r, pcl_z, pcl_var\n")
    for pos_id in range(var_num):
        pcl_r.append(proc_pcl_r[pos_id])
        pcl_z.append(proc_pcl_z[pos_id])
        pcl_var.append(proc_pcl_var[pos_id])
        # output data, pcl_r, pcl_z, pcl_var
        csv_file.write("%f, %f, %f,\n" %  (pcl_r[pos_id], pcl_z[pos_id], pcl_var[pos_id]))
    csv_file.close()

    # plot data
    fig = plt.figure()
    plot1 = fig.subplots(1, 1)
    line1, = plot1.plot(pcl_r, pcl_var)
    plt.show()
