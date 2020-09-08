import math
import h5py as py

def plate_with_hole_stress(x, y, r, load):
    rou = math.sqrt(x*x + y*y)
    fai = math.atan(y/x)
    r2 = r * r
    rou2 = rou * rou
    cos_2fai = math.cos(2.0*fai)
    sin_2fai = math.sin(2.0*fai)
    s_rou = 0.5 * load * (1.0 - r2/rou2 + cos_2fai * (1.0 - r2 / rou2) * (1.0 - 3.0 * r2 / rou2))
    s_fai = 0.5 * load * (1.0 + r2/rou2 - cos_2fai * (1.0 + 3.0 * r2*r2 / (rou2*rou2)))
    s_roufai = -0.5 * load * sin_2fai * (1.0 - r2/rou2) * (1.0 + 3.0 * r2 / rou2)
    sin_fai = math.sin(fai)
    cos_fai = math.cos(fai)
    sin2_fai = sin_fai * sin_fai
    cos2_fai = cos_fai * cos_fai
    sincos_fai = sin_fai * cos_fai
    s11 = s_rou * cos2_fai + s_fai * sin2_fai - 2.0 * s_roufai * sincos_fai
    s22 = s_rou * sin2_fai + s_fai * cos2_fai + 2.0 * s_roufai * sincos_fai
    s12 = (s_rou - s_fai) * sincos_fai + s_roufai * (cos2_fai - sin2_fai)
    return s11, s22, s12

hole_r = 1.0
surface_load = -10.0

res_file = open("plate_with_hole_res.csv", "w")
res_file.write("x, y, cal_s11, acc_s11, s11_error, cal_s22, acc_s22, s22_error, cal_s12, acc_s12, s12_error,\n")

hdf5_file = py.File("../Build/Tests/t2d_me_s_plate_with_hole.h5", "r")
th_grp = hdf5_file['TimeHistory']['load1']
output_num = th_grp.attrs['output_num']
# last frame
frame_grp = th_grp['frame_%d' % (output_num-1)]
pcl_data_grp = frame_grp['ParticleData']
pcl_num = pcl_data_grp.attrs['pcl_num']
pcl_dset = pcl_data_grp['field']

for pcl_id in range(pcl_num):
    pcl_fld = pcl_dset[pcl_id]
    pcl_x = pcl_fld['x']
    pcl_y = pcl_fld['y']
    if pcl_x <= 5.0 and pcl_y <= 5.0:
        cal_s11 = pcl_fld['s11']
        cal_s22 = pcl_fld['s22']
        cal_s12 = pcl_fld['s12']
        acc_s11, acc_s22, acc_s12 = plate_with_hole_stress(pcl_x, pcl_y, hole_r, surface_load)
        if abs(acc_s11) != 0.0:
            s11_error = abs(cal_s11 - acc_s11) / abs(acc_s11) * 100.0
        else:
            s11_error = 100.0
        if abs(acc_s22) != 0.0:
            s22_error = abs(cal_s22 - acc_s22) / abs(acc_s22) * 100.0
        else:
            s22_error = 100.0
        if abs(acc_s12) != 0.0:
            s12_error = abs(cal_s12 - acc_s12) / abs(acc_s12) * 100.0
        else:
            s12_error = 100.0
        res_file.write("%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f,\n" % \
            (pcl_x, pcl_y, cal_s11, acc_s11, s11_error, cal_s22, acc_s22, s22_error, cal_s12, acc_s12, s12_error))

hdf5_file.close()
res_file.close()
