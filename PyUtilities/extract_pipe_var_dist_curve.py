import math
import sys
import h5py as py
import matplotlib.pyplot as plt
from multiprocessing import Process, Array

path_name = "../Build/TestsParallel/"
sim_name = "t2d_chm_mt_pipe_embedment_geo_rou_r05"
var_name = "p"
th_name = "geostatic"
rb_name = "RigidCircle"
th_id = 101
proc_num = 2

# (x0, y0, dx, dy)
# absolute rb_position
var_pos = [
    (0.1, -0.1, 0.2, 0.2),
    (0.1, -0.3, 0.2, 0.2),
    (0.1, -0.5, 0.2, 0.2),
    (0.1, -0.75, 0.2, 0.2),
    (0.1, -1.0, 0.2, 0.2),
    (0.1, -1.25, 0.2, 0.2),
    (0.1, -1.5, 0.2, 0.2),
    (0.1, -1.7, 0.2, 0.2),
    (0.1, -2.0, 0.2, 0.2),
    (0.1, -2.5, 0.2, 0.2),
    (0.1, -3.0, 0.2, 0.2),
    (0.1, -4.0, 0.2, 0.2),
    (0.1, -6.0, 0.2, 0.2),
    (0.1, -8.0, 0.2, 0.2),
]

###########################################################################################
###########################################################################################
def get_avg_var(pcl_xs, pcl_ys, pcl_vols, pcl_vars, pcl_num,  \
                model_var_x, model_var_y, model_var_dx, model_var_dy):
    avg_var = 0.0
    eff_pcl_vol = 0.0
    eff_pcl_num = 0
    model_var_xmin = model_var_x - model_var_dx
    model_var_xmax = model_var_x + model_var_dx
    model_var_ymin = model_var_y - model_var_dy
    model_var_ymax = model_var_y + model_var_dy
    for pcl_id in range(pcl_num):
        p_x = pcl_xs[pcl_id]
        p_y = pcl_ys[pcl_id]
        if p_x >= model_var_xmin and p_x <= model_var_xmax and \
            p_y >= model_var_ymin and p_y <= model_var_ymax:
            avg_var += pcl_vars[pcl_id] * pcl_vols[pcl_id]
            eff_pcl_vol += pcl_vols[pcl_id]
            eff_pcl_num += 1
    if (eff_pcl_vol != 0.0):
        avg_var /= eff_pcl_vol
    # average pore pressure, effective pcl vol, effective pcl num
    return (avg_var, eff_pcl_vol, eff_pcl_num)

class AvgVarParam:
    def __init__(self, path_name, sim_name, th_name, th_id, \
                 rb_name, var_name, rb_ini_y, var_pos):
        self.path_name = path_name
        self.sim_name = sim_name
        self.th_name = th_name
        self.th_id = th_id
        self.rb_name = rb_name
        self.var_name = var_name
        self.rb_ini_y = rb_ini_y
        self.var_pos = var_pos

def proc_get_th_avg_var(pos_id_range, param, pcl_x, pcl_y, pcl_var):
    hdf5_file = py.File(param.path_name + param.sim_name + ".h5", "r")
    th_grp = hdf5_file['TimeHistory'][param.th_name]
    frame_grp = th_grp['frame_%d' % param.th_id]
    rb_grp = frame_grp[param.rb_name]
    rb_dy = rb_grp.attrs['y'] - param.rb_ini_y
    # cal avg pore pressure for each time histories
    for pos_id in range(pos_id_range[0], pos_id_range[1]):
        #
        avg_x = param.var_pos[pos_id][0]
        avg_y = param.var_pos[pos_id][1] + rb_dy
        # pcl x, y pos results
        pcl_x[pos_id] = avg_x
        pcl_y[pos_id] = avg_y
        pcl_dset = frame_grp['ParticleData']['field']
        pcl_xs = pcl_dset[:]['x']
        pcl_ys = pcl_dset[:]['y']
        pcl_vols = pcl_dset[:]['vol']
        pcl_vars = pcl_dset[:][param.var_name]
        pcl_num = len(pcl_dset)
        avg_var, eff_pcl_vol, eff_pcl_num = get_avg_var( \
            pcl_xs, pcl_ys, pcl_vols, pcl_vars, pcl_num, \
            avg_x, avg_y, param.var_pos[pos_id][2], param.var_pos[pos_id][3])
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
    rb_ini_y = rb_grp.attrs['y']
    hdf5_file.close()
    
    var_num = len(var_pos)
    proc_pcl_x = Array('d', [0.0] * var_num)
    proc_pcl_y = Array('d', [0.0] * var_num)
    proc_pcl_var = Array('d', [0.0] * var_num)
    procs = []
    proc_param = AvgVarParam(path_name, sim_name, th_name, th_id, rb_name, var_name, rb_ini_y, var_pos)
    for proc_id in range(proc_num):
        p = Process(target = proc_get_th_avg_var, \
                    args = ([int(proc_id*var_num/proc_num), int((proc_id+1)*var_num/proc_num)], \
                             proc_param, proc_pcl_x, proc_pcl_y, proc_pcl_var))
        procs.append(p)
        p.start()
    
    print("All procs started.")
    
    for proc_id in range(proc_num):
        procs[proc_id].join()
        print("proc %d completed." % proc_id)
    
    print("All procs completed.")

    pcl_x = []
    pcl_y = []
    pcl_var = []
    csv_file = open(path_name + sim_name + "_" + var_name + ".csv", "w")
    csv_file.write("pcl_x, pcl_y, pcl_var\n")
    for pos_id in range(var_num):
        pcl_x.append(proc_pcl_x[pos_id])
        pcl_y.append(proc_pcl_y[pos_id])
        pcl_var.append(proc_pcl_var[pos_id])
        # output data, pcl_x, pcl_y, pcl_var
        csv_file.write("%f, %f, %f,\n" %  (pcl_x[pos_id], pcl_y[pos_id], pcl_var[pos_id]))
    csv_file.close()

    # plot data
    fig = plt.figure()
    plot1 = fig.subplots(1, 1)
    line1, = plot1.plot(pcl_y, pcl_var)
    plt.show()
