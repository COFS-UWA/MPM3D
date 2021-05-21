import h5py as py
import matplotlib.pyplot as plt

# Numerical result
hdf5_file = py.File("../Build/TestsParallel/t3d_me_mt_piezofoundation_geo.h5", "r")

th_grp = hdf5_file['TimeHistory']['geostatic']
th_num = th_grp.attrs['output_num']

sim_times = []
f_ubs = []
e_kins = []
is_init = False
for th_id in range(th_num):
    th_gp = th_grp['frame_%d' % th_id]
    sim_t = th_gp.attrs['total_time']
    f_ub = th_gp.attrs['f_ub_ratio']
    e_kin = th_gp.attrs['e_kin_ratio']
    sim_times.append(sim_t)
    f_ubs.append(f_ub)
    e_kins.append(e_kin)

hdf5_file.close()

data_file = open("../Build/TestsParallel/t3d_me_mt_piezofoundation_geo_ratio.csv", "w")
for i in range(len(sim_times)):
    data_file.write("%f, %f, %f\n" % (sim_times[i], f_ubs[i], e_kins[i]))
data_file.close()

print(e_kins)

fig = plt.figure()
plot1 = fig.subplots(1, 1)
line1, = plot1.plot(sim_times, e_kins)

sim_time_range = [min(sim_times), max(sim_times)]
plt.xlim(sim_time_range)
#plt.ylim(y_range)

plt.show()
