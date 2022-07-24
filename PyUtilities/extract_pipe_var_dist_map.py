import math
import sys
import numpy as np
import h5py as py
import matplotlib.pyplot as plt
import csv
from multiprocessing import Process, Array

var_pos_csv = "var_pos"

path_name = "../Build/TestsParallel/"
sim_name = "t2d_chm_mt_pipe_embedment_geo_rou_r05"
var_name = "s22"
th_name = "geostatic"
th_id = 101
proc_num = 2

# (x0, y0, dx, dy)
# absolute rb_position
var_dx = 0.2
var_dy = 0.2

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
                 var_name, var_pos, var_dx, var_dy):
        self.path_name = path_name
        self.sim_name = sim_name
        self.th_name = th_name
        self.th_id = th_id
        self.var_name = var_name
        self.var_pos = var_pos
        self.var_dx = var_dx
        self.var_dy = var_dy

def proc_get_th_avg_var(pos_id_range, param, pcl_x, pcl_y, pcl_var):
    hdf5_file = py.File(param.path_name + param.sim_name + ".h5", "r")
    th_grp = hdf5_file['TimeHistory'][param.th_name]
    frame_grp = th_grp['frame_%d' % param.th_id]
    # cal avg pore pressure for each time histories
    for pos_id in range(pos_id_range[0], pos_id_range[1]):
        #
        avg_x = param.var_pos[pos_id][0]
        avg_y = param.var_pos[pos_id][1]
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
            avg_x, avg_y, param.var_dx, param.var_dy)
        pcl_var[pos_id] = avg_var
        # print msg
        if pos_id_range[0] == 0:
            print("processed %d%%, %d pcls nearby." % (int((pos_id+1)/(pos_id_range[1]-pos_id_range[0])*100.0), eff_pcl_num))
    hdf5_file.close()

if __name__ == "__main__":
    out_file = open(path_name + sim_name + "_" + var_name + ".csv", "w")
    csv_file = open(var_pos_csv + ".csv", "r")
    var_pos_reader = csv.reader(csv_file)
    
    # read x coordinates
    x_var_pos_str = next(var_pos_reader)
    x_var_num = len(x_var_pos_str)
    x_var_pos = np.zeros(x_var_num)
    for i in range(1, x_var_num):
        x_var_pos[i] = float(x_var_pos_str[i])
        out_file.write("%f, " % x_var_pos[i])
    out_file.write("\n")

    # read y coordinates
    y_var_pos_str = next(var_pos_reader)
    y_var_num = len(y_var_pos_str)
    y_var_pos = np.zeros(y_var_num)
    for i in range(1, y_var_num):
        y_var_pos[i] = float(y_var_pos_str[i])
        out_file.write("%f, " % y_var_pos[i])
    out_file.write("\n")

    # read coordinates
    var_pos = []
    var_pos_line_num = 0
    for row in var_pos_reader:
        row2 = next(var_pos_reader)
        row_coords = []
        for i in range(1, len(row)):
            if row[i] != '' and row2[i] != '':
                row_coords.append((float(row[i]), float(row2[i])))
        var_pos.append(row_coords)
        var_pos_line_num += 1
    
    # extract variables
    for row_id in range(len(var_pos)):
        var_row_pos = var_pos[row_id]
        
        var_num = len(var_row_pos)
        proc_pcl_x = Array('d', [ 0.0 ] * var_num)
        proc_pcl_y = Array('d', [ 0.0 ] * var_num)
        proc_pcl_var = Array('d', [ 0.0 ] * var_num)
        procs = []
        proc_param = AvgVarParam(path_name, sim_name, th_name, th_id, \
                                 var_name, var_row_pos, var_dx, var_dy)
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
        print("Complete row %d" % row_id)

        for pos_id in range(var_num):
            out_file.write("%f," % proc_pcl_x[pos_id])
        out_file.write("\n")
        for pos_id in range(var_num):
            out_file.write("%f," % proc_pcl_y[pos_id])
        out_file.write("\n")
        for pos_id in range(var_num):
            out_file.write("%f," % proc_pcl_var[pos_id])
        out_file.write("\n")
    
    out_file.close()
    csv_file.close()
