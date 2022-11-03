import numpy as np
import matplotlib.pyplot as plt
import h5py as py

from OneDConsolidation import OneDConsolidation

fig = plt.figure()
plot1 = fig.subplots(1, 1)
plot1.set_xlabel("time")
plot1.set_ylabel("displacement")

# analytical solution
u0 = 1.0
H = 1.0
E = 1000.0
niu = 0.0 # possion ratio
kv = 1.0e-4
miu = 1.0 # dynamic viscosity

Es = (1 - niu) / (1 + niu) / (1 - 2.0*niu) * E # Es = (1-v) / (1 + v) / (1-2v) * E
Cv = kv * Es / miu
con_res = OneDConsolidation(Cv, Es, u0, H)
time = 10.0 # time of consolidation
data_num = 100
t_list = np.zeros(data_num + 2)
u_list = np.zeros(data_num + 2)
t_list[0] = 0.0
u_list[0] = 0.0
t_list[1] = 0.0 # time for equilibrium
u_list[1] = u_list[0]
for i in range(data_num):
    t_list[i + 2] = time * float(i) / float(data_num)
    u_list[i + 2] = con_res.calSettlement(t_list[i + 2])
    t_list[i + 2] += t_list[1]

plot1.plot(t_list, u_list, 'r--')
with open("consolidation_disp_ana_SE.csv", "w") as out_file:
    for i in range(len(t_list)):
        out_file.write("%f, %f\n" % (t_list[i], u_list[i]))

plt.show()