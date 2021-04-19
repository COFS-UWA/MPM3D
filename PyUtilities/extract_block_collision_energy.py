import h5py as py
import numpy as np
import matplotlib.pyplot as plt

hdf5_file = py.File("../Build/TestsParallel/t2d_me_mt_block_collision.h5", "r")

th_grp = hdf5_file['TimeHistory']['collision']
th_num = th_grp.attrs['output_num']
pcl_num = th_grp['frame_0']['ParticleData'].attrs['pcl_num']


E = 10000.0
niu = 0.0
lmbd2 = E * niu / ((1.0 + niu) * (1.0 - 2.0 * niu))
G = E / (1.0 + niu)
lmbd1 = lmbd2 + G
sim_time = []
db_kinetic_energy = []
db_strain_energy = []
rb_kinetic_energy = []
db_energy = []
kinetic_energy = []
total_energy = []
prev_p_s = np.zeros([pcl_num, 3])
prev_p_e = np.zeros([pcl_num, 3])
for th_id in range(th_num):
    frame_grp = th_grp['frame_%d' % th_id]
    sim_time.append(frame_grp.attrs['total_time'])

    db_dset = frame_grp['ParticleData']['field']
    db_kin = 0.0
    db_se = 0.0
    for p_id in range(pcl_num):
        p_var = db_dset[p_id]
        p_vx = p_var['vx']
        p_vy = p_var['vy']
        p_m = p_var['m']
        p_vol = p_m / p_var['density']
        db_kin += 0.5 * p_m * (p_vx * p_vx + p_vy * p_vy) 
        p_s11 = p_var['s11']
        p_s22 = p_var['s22']
        p_s12 = p_var['s12']
        p_e11 = p_var['e11']
        p_e22 = p_var['e22']
        p_e12 = p_var['e12']
        # db_se += p_vol * 0.5 * \
                # ((prev_p_s[p_id][0] + p_s11) * (p_e11 - prev_p_e[p_id][0]) \
               # + (prev_p_s[p_id][1] + p_s22) * (p_e22 - prev_p_e[p_id][1])
               # + (prev_p_s[p_id][2] + p_s12) * (p_e12 - prev_p_e[p_id][2]))
        db_se += p_vol * (0.5 * lmbd1 * (p_e11 * p_e11 + p_e22 * p_e22) \
               + lmbd2 * p_e11 * p_e22 + 0.5 * G * p_e12 * p_e12)
        prev_p_s[p_id][0] = p_s11
        prev_p_s[p_id][1] = p_s22
        prev_p_s[p_id][2] = p_s12
        prev_p_e[p_id][0] = p_e11
        prev_p_e[p_id][1] = p_e22
        prev_p_e[p_id][2] = p_e12
    db_kinetic_energy.append(db_kin)
    db_strain_energy.append(db_se)
    
    rb_grp = frame_grp['RigidRect']
    rb_vx = rb_grp.attrs['vx']
    rb_vy = rb_grp.attrs['vy']
    rb_vang = rb_grp.attrs['v_angle']
    rb_den = rb_grp.attrs['density']
    rb_hx = rb_grp.attrs['hx']
    rb_hy = rb_grp.attrs['hy']
    rb_m = rb_hx * rb_hy * rb_den
    rb_moi = rb_m * (rb_hx * rb_hx + rb_hy * rb_hy) / 12.0
    rb_kin = 0.5 * rb_m * (rb_vx * rb_vx + rb_vy * rb_vy) \
           + 0.5 * rb_moi * rb_vang * rb_vang
    rb_kinetic_energy.append(rb_kin)
    
    db_energy.append(db_kinetic_energy[th_id] + db_strain_energy[th_id])
    kinetic_energy.append(db_kinetic_energy[th_id] + rb_kinetic_energy[th_id])
    total_energy.append(kinetic_energy[th_id] + db_strain_energy[th_id])

hdf5_file.close()

data_file = open("../Build/TestsParallel/t2d_me_mt_block_collision_energy.csv", "w")
data_file.write("time, db_kin, db_seg, db_total, rb_kin, kin, total\n")
for i in range(th_num):
    data_file.write("%f, %f, %f, %f, %f, %f, %f\n" % \
        (sim_time[i], db_kinetic_energy[i], db_strain_energy[i], db_energy[i], \
         rb_kinetic_energy[i], kinetic_energy[i], total_energy[i]))
data_file.close()

fig = plt.figure()
plot1 = fig.subplots(1, 1)
line1, = plot1.plot(sim_time, total_energy)
line2, = plot1.plot(sim_time, rb_kinetic_energy)
line3, = plot1.plot(sim_time, db_kinetic_energy)
line4, = plot1.plot(sim_time, db_strain_energy)

x_range = [min(sim_time), max(sim_time)]
plt.xlim(x_range)

plt.legend(handles=[ line1, line2, line3, line4 ], \
    labels=[ 'Total', 'rb_kinetic', 'db_kinetic', 'db_strain' ])
plt.show()
