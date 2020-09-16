import numpy as np
import matplotlib.pyplot as plt
import h5py as py

pcl_id = 54337
frame_num = 7

hdf5_file = py.File("../Build/Tests/t2d_chm_s_pipe_conference_restart3.h5", "r")
th_grp = hdf5_file['TimeHistory']['penetration']

frame_ids = []
pcl_pcs = []

for frame_id in range(frame_num):
    frame_grp = th_grp['frame_%d' % frame_id]
    pcl_grp = frame_grp['ParticleData']
    pcl_dset = pcl_grp['field']
    mcc_grp = frame_grp['MaterialModel']
    mcc_dset = mcc_grp['ModifiedCamClay']
    
    pcl_data = pcl_dset[pcl_id]
    mcc_data = mcc_dset[pcl_id]
    
    frame_ids.append(float(frame_id))
    pcl_pcs.append(float(mcc_data['s22']))

hdf5_file.close()

fig = plt.figure()
plot1 = fig.subplots(1, 1)
plot1.set_xlabel("frame")
plot1.set_ylabel("pc")
line1, = plot1.plot(frame_ids, pcl_pcs)
plt.show()
