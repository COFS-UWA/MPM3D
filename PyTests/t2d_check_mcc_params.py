import numpy as np
import matplotlib.pyplot as plt
import h5py as py

t_id = 424

hdf5_file = py.File("../Build/Tests/t2d_chm_s_pipe_conference_restart2.h5", "r")
th_grp = hdf5_file['TimeHistory']['penetration']

frame_grp = th_grp['frame_%d' % t_id]
mcc_grp = frame_grp['MaterialModel']
mcc_dset = mcc_grp['ModifiedCamClay']
pcl_num = len(mcc_dset)

# pc_min_pcl_id = 0
# pc_min = mcc_dset[0]['pc']
# pc_max_pcl_id = 0
# pc_max = pc_min
# for pcl_id in range(1, pcl_num):
    # mcc_data = mcc_dset[pcl_id]
    # pc = mcc_data['pc']
    # if (pc_min > pc):
        # pc_min = pc
        # pc_min_pcl_id = pcl_id
    # if (pc_max < pc):
        # pc_max = pc
        # pc_max_pcl_id = pcl_id


# print("frame %d\nmin %.2f at pcl %d\nmax %.2f at pcl %d\n" \
    # % (t_id, pc_min, pc_min_pcl_id, pc_max, pc_max_pcl_id))

pcl_id = 54337
pcl_pc = mcc_dset[pcl_id]['pc']
print("frame %d pc %.2f at pcl %d\n" % (t_id, pcl_pc, pcl_id))

hdf5_file.close()
