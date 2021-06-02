import math
import matplotlib.pyplot as plt

phi = 30.0
p = 10.0
lode = 30.0

phi = math.radians(phi)
sin_phi = math.sin(phi)
sqrt3 = math.sqrt(3.0)

lode_angs = []
q_mcs = []
q_mns = []
phi_equs = []
lode_ang_degs = []
intv_num = 100
lode_intv = math.pi / (3.0 * float(intv_num))
for i in range(-intv_num + 1, intv_num + 1, 1):
    lode_angs.append(i * lode_intv)
    lode_ang_degs.append(math.degrees(i * lode_intv))

p_sin_phi = -6.0 * p * sin_phi
A = (1.0 - sin_phi*sin_phi) / (9.0 - sin_phi*sin_phi)
C = (A - 1.0/3.0) * p
D = (1.0 - 9.0*A) * p*p*p
for la in lode_angs:
    sin_la = math.sin(la)
    cos_la = math.cos(la)
    q_mc = 0.0
    if la >= 0.0:
        q_mc = p_sin_phi / (sin_phi*(cos_la-sqrt3*sin_la)-sqrt3*(sqrt3*cos_la+sin_la))
    else:    
        q_mc = p_sin_phi / (sin_phi*(cos_la+sqrt3*sin_la)-sqrt3*(sqrt3*cos_la-sin_la))
    q_mcs.append(q_mc)
    B = 2.0 / 27.0 * math.cos(3.0*la)
    qn = q_mc
    qn_1 = 0.0
    while (True):
        qn_1 = qn - (B*qn*qn*qn + C*qn*qn + D) / (3.0*B*qn*qn + 2.0*C*qn)
        if (qn_1 - qn) / qn < 0.001:
            break
        qn = qn_1
    q_mns.append(qn_1)
    phi_equ = 0.0
    if la >= 0.0:
        phi_equ = sqrt3*qn_1*(sqrt3*cos_la+sin_la)/(6.0*p+qn_1*(cos_la-sqrt3*sin_la))
    else:
        phi_equ = sqrt3*qn_1*(sqrt3*cos_la-sin_la)/(6.0*p+qn_1*(cos_la+sqrt3*sin_la))
    phi_equs.append(math.degrees(math.asin(phi_equ))) # asin
    
print(phi_equs)

fig = plt.figure()

fig.add_subplot(1, 2, 1, projection='polar')
# step4: 绘制极坐标图
plt.polar(lode_angs, q_mcs)
plt.polar(lode_angs, q_mns)

fig.add_subplot(1, 2, 2)
plt.plot(lode_ang_degs, phi_equs)

plt.show()