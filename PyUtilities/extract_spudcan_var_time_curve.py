import math
import sys
import h5py as py
import matplotlib.pyplot as plt
from multiprocessing import Process, Array

path_name = "E:/"
#path_name = "../Build/TestsParallel/"
sim_name = "t3d_chm_mt_spudcan_cy_ucav-130_V2000"
var_name = "p"
model_var_r0 = 0.1
model_var_z0 = -0.2
model_var_dr = 0.2
model_var_dz = 0.2
proc_num = 4

# time history name
th_name = "penetration"
rb_name = "RigidObjectByT3DMesh"

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
    def __init__(self, path_name, sim_name, th_name, rb_name, var_name, \
                 rb_ini_z, model_var_r0, model_var_z0, \
                 model_var_dr, model_var_dz, \
                 model_cen_x = 0.0, model_cen_y = 0.0):
        self.path_name = path_name
        self.sim_name = sim_name
        self.th_name = th_name
        self.rb_name = rb_name
        self.var_name = var_name
        self.rb_ini_z = rb_ini_z
        self.model_var_r0 = model_var_r0
        self.model_var_z0 = model_var_z0
        self.model_var_dr = model_var_dr
        self.model_var_dz = model_var_dz
        self.model_cen_x = model_cen_x
        self.model_cen_y = model_cen_y

def proc_get_th_avg_var(th_id_range, param, \
    rb_z, pcl_r, pcl_z, pcl_var):
    #
    hdf5_file = py.File(param.path_name + param.sim_name + ".h5", "r")
    th_grp = hdf5_file['TimeHistory'][param.th_name]
    th_num = th_grp.attrs['output_num']
    # cal avg pore pressure for each time histories
    for th_id in range(th_id_range[0], th_id_range[1]):
        frame_grp = th_grp['frame_%d' % th_id]
        rb_grp = frame_grp[param.rb_name]
        rb_dz = rb_grp.attrs['z'] - param.rb_ini_z
        rb_z[th_id] = -rb_dz
        #
        avg_r = param.model_var_r0
        avg_z = param.model_var_z0 + rb_dz
        pcl_r[th_id] = avg_r
        pcl_z[th_id] = avg_z
        pcl_dset = frame_grp['ParticleData']['field']
        pcl_xs = pcl_dset[:]['x']
        pcl_ys = pcl_dset[:]['y']
        pcl_zs = pcl_dset[:]['z']
        pcl_vols = pcl_dset[:]['vol']
        pcl_vars = pcl_dset[:][param.var_name]
        pcl_num = len(pcl_dset)
        avg_var, eff_pcl_vol, eff_pcl_num = get_avg_var( \
            pcl_xs, pcl_ys, pcl_zs, pcl_vols, pcl_vars, pcl_num, \
            avg_r, avg_z, param.model_var_dr, param.model_var_dz, \
            param.model_cen_x, param.model_cen_y)
        pcl_var[th_id] = avg_var
        # print msg
        if th_id_range[0] == 0:
            print("processed %d%%, %d pcls nearby." % (int((th_id+1)/(th_id_range[1]-th_id_range[0])*100.0), eff_pcl_num))
    # complete
    hdf5_file.close()

if __name__ == "__main__":
    hdf5_file = py.File(path_name + sim_name + ".h5", "r")
    th_grp = hdf5_file['TimeHistory'][th_name]
    th_num = th_grp.attrs['output_num']
    frame_grp = th_grp['frame_0']
    rb_grp = frame_grp[rb_name]
    rb_ini_z = rb_grp.attrs['z']
    hdf5_file.close()
    
    proc_rb_z = Array('d', [0.0]*th_num)
    proc_pcl_r = Array('d', [0.0]*th_num)
    proc_pcl_z = Array('d', [0.0]*th_num)
    proc_pcl_var = Array('d', [0.0]*th_num)
    procs = []
    proc_param = AvgVarParam(path_name, sim_name, th_name, rb_name, var_name, \
                             rb_ini_z, model_var_r0, model_var_z0, model_var_dr,model_var_dz)
    for proc_id in range(proc_num):
        p = Process(target = proc_get_th_avg_var, \
                    args = ([int(proc_id * th_num / proc_num), int((proc_id+1) * th_num / proc_num)], \
                             proc_param, proc_rb_z, proc_pcl_r, proc_pcl_z, proc_pcl_var))
        procs.append(p)
        p.start()
    
    print("All procs started.")
    
    for proc_id in range(proc_num):
        procs[proc_id].join()
        print("proc %d completed." % proc_id)
    
    print("All procs completed.")

    rb_z = []
    pcl_r = []
    pcl_z = []
    pcl_var = []
    csv_file = open(path_name + sim_name + "_" + var_name + ".csv", "w")
    csv_file.write("rb_z, pcl_r, pcl_z, pcl_var\n")
    for th_id in range(th_num):
        rb_z.append(proc_rb_z[th_id])
        pcl_r.append(proc_pcl_r[th_id])
        pcl_z.append(proc_pcl_z[th_id])
        pcl_var.append(proc_pcl_var[th_id])
        # output data
        # rb_z, pcl_r, pcl_z, pcl_var
        csv_file.write("%f, %f, %f, %f,\n" %  (rb_z[th_id], pcl_r[th_id], pcl_z[th_id], pcl_var[th_id]))
    csv_file.close()

    # plot data
    fig = plt.figure()
    plot1 = fig.subplots(1, 1)
    line1, = plot1.plot(pcl_var, rb_z)
    plt.ylim(rb_z[0], rb_z[-1])
    plt.show()
