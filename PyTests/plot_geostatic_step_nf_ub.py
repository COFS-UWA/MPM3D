import numpy as np
import matplotlib.pyplot as plt
import h5py as py

fig = plt.figure()
plot1 = fig.subplots(1, 1)
plot1.set_xlabel("time")
plot1.set_ylabel("kin_energy")

out_time = []
max_nf_ub = []
nf_ratio = []
kin_energy = []

# numerical solution
hdf5_file = py.File("..\\Build\\Tests\\t2d_me_s_geostatic.h5", "r")
th_grp = hdf5_file['TimeHistory']['geostatic']
output_num = th_grp.attrs['output_num']
for t_id in range(output_num):
    # frame
    frame_grp = th_grp['frame_%d' % t_id]
    frame_time = frame_grp.attrs['total_time']
    out_time.append(frame_time)
    mfb = frame_grp.attrs['max_unbalanced_nodal_force']
    max_nf_ub.append(mfb)
    nfr = frame_grp.attrs['unbalanced_nodal_force_ratio']
    nf_ratio.append(nfr)
    ke = frame_grp.attrs['kinetic_energy']
    kin_energy.append(ke)

hdf5_file.close()

#line1, = plot1.plot(out_time, max_nf_ub)
#line1, = plot1.plot(out_time, nf_ratio)
line1, = plot1.plot(out_time, kin_energy)
plt.show()
