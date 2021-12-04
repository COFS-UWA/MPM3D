import h5py as py
import matplotlib.pyplot as plt

# Numerical result
hdf5_file = py.File("../Build/TestsParallel/t2d_chm_mt_block_sliding.h5", "r")

th_grp = hdf5_file['TimeHistory']['slide']
th_num = th_grp.attrs['output_num']

is_init = False
ini_x = 0.0
rb_x = []
rb_fsx_cont = []
rb_fsy_cont = []
rb_ffx_cont = []
rb_ffy_cont = []
rb_fx_cont = []
rb_fy_cont = []
for th_id in range(th_num):
    rb_grp = th_grp['frame_%d' % th_id]['RigidRect']
    cen_x = rb_grp.attrs['x']
    rb_fsx = rb_grp.attrs['fsx_cont']
    rb_fsy = rb_grp.attrs['fsy_cont']
    rb_ffx = rb_grp.attrs['ffx_cont']
    rb_ffy = rb_grp.attrs['ffy_cont']
    rb_fx = rb_grp.attrs['fx_contact']
    rb_fy = rb_grp.attrs['fy_contact']
    if not is_init:
        is_init = True
        ini_x = cen_x
    else:
        rb_x.append(ini_x - cen_x)
        rb_fsx_cont.append(rb_fsx)
        rb_fsy_cont.append(rb_fsy)
        rb_ffx_cont.append(rb_ffx)
        rb_ffy_cont.append(rb_ffy)
        rb_fx_cont.append(rb_fx)
        rb_fy_cont.append(rb_fy)

hdf5_file.close()

data_file = open("../Build/TestsParallel/t2d_chm_mt_block_sliding_rf.csv", "w")
data_file.write("x_disp, fsx_cont, fsy_cont, ffx_cont, ffy_cont, fx_cont, fy_cont,\n")
for i in range(len(rb_x)):
    data_file.write("%f, %f, %f, %f, %f, %f, %f,\n" % (rb_x[i], \
        rb_fsx_cont[i], rb_fsy_cont[i], \
        rb_ffx_cont[i], rb_ffy_cont[i], \
        rb_fx_cont[i], rb_fy_cont[i]))
data_file.close()

fig = plt.figure()
plot1 = fig.subplots(1, 1)
line1, = plot1.plot(rb_x, rb_fsy_cont)

plt.show()
