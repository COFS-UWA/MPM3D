import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import math as mh
import matplotlib.pyplot as plt
import csv

path_name = "../Build/TestsParallel/"
var_pos_csv = "t2d_chm_mt_pipe_embedment_geo_rou_r05_s22"

# a value indicating meaningless number
padding_value = 1.0e50

def find_index(x, x_array):
    max_index = len(x_array) - 1
    if x < x_array[0]:
        return 0
    if x > x_array[max_index]:
        return max_index
    mid_index = mh.floor(max_index/2)
    min_index = 0
    while min_index != mid_index:
        if x < x_array[mid_index]:
            max_index = mid_index
        else:
            min_index = mid_index
        mid_index = mh.floor((min_index + max_index)/2)
    if x - x_array[min_index] < x_array[mid_index] - x:
        return min_index
    return mid_index

# 设置图例字号
mpl.rcParams['legend.fontsize'] = 10

# 方式2：设置三维图形模式
fig = plt.figure()
ax = fig.gca(projection='3d')

csvfile = open(path_name + var_pos_csv + ".csv", newline='')
var_pos_reader = csv.reader(csvfile)

# read x coordinates
x_var_pos_str = next(var_pos_reader)
x_var_num = 0
for pos in x_var_pos_str:
    if pos != " ":
        x_var_num += 1
x_var_pos = np.zeros(x_var_num)
for i in range(x_var_num):
    x_var_pos[i] = float(x_var_pos_str[i])
x_var_pos.sort()

# read y coordinates
y_var_pos_str = next(var_pos_reader)
y_var_num = 0
for pos in y_var_pos_str:
    if pos != " ":
        y_var_num += 1
y_var_pos = np.zeros(y_var_num)
for i in range(y_var_num):
    y_var_pos[i] = float(y_var_pos_str[i])
y_var_pos.sort()

# read variable values
hor_var_value = np.zeros([y_var_num, x_var_num])
ver_var_value = np.zeros([x_var_num, y_var_num])
for x_id in range(x_var_num):
    for y_id in range(y_var_num):
        hor_var_value[y_id][x_id] = padding_value
        ver_var_value[x_id][y_id] = padding_value

for row_x in var_pos_reader:
    row_y = next(var_pos_reader)
    row_var = next(var_pos_reader)
    for i in range(len(row_x)):
        if row_x[i] != '' and row_y[i] != '' and row_var[i] != '':
            x_pos = float(row_x[i])
            y_pos = float(row_y[i])
            var = float(row_var[i])
            x_index = find_index(x_pos, x_var_pos)
            y_index = find_index(y_pos, y_var_pos)
            hor_var_value[y_index][x_index] = var
            ver_var_value[x_index][y_index] = var

csvfile.close()

print(y_var_num)

for y_id in range(y_var_num):
    var_x = []
    var_y = []
    var_value = []
    for x_id in range(x_var_num):
        if hor_var_value[y_id][x_id] != padding_value:
            var_x.append(x_var_pos[x_id])
            var_y.append(y_var_pos[y_id])
            var_value.append(hor_var_value[y_id][x_id])
    ax.plot(var_x, var_y, var_value, color='blue')

for x_id in range(x_var_num):
    var_x = []
    var_y = []
    var_value = []
    for y_id in range(y_var_num):
        if ver_var_value[x_id][y_id] != padding_value:
            var_x.append(x_var_pos[x_id])
            var_y.append(y_var_pos[y_id])
            var_value.append(ver_var_value[x_id][y_id])
    ax.plot(var_x, var_y, var_value, color='blue')

# 显示图例
#ax.legend()

# 显示图形
plt.show()
