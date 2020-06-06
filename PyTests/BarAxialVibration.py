import math
import numpy as np
import matplotlib.pyplot as plt

class BarAxialVibration:
    """
    1 dimensional bar vibration.
    """
    def __init__(self, H, p0, bf, E, density, error_ratio = 1.0e-3):
        self.H = H
        self.p0 = p0
        self.bf = bf
        self.E = E
        self.density = density
        self.error_ratio = error_ratio
        self.alpha = math.sqrt(self.E/self.density)
    
    def displacement(self, y, t):
        u = 0.0
        sign = 1.0
        i = 0
        max_iter_num = 100
        while (True):
            sign *= -1.0
            i += 1
            u_n = 8.0*self.H * (2.0*math.pi*self.p0*i*sign - 2.0*self.density*self.bf*self.H - math.pi*self.p0*sign) / ((4.0*i*i - 4.0*i + 1.0) * (2.0*i - 1.0) * math.pi*math.pi*math.pi * self.E)
            tmp = u_n * math.cos(self.alpha*(2.0*i-1.0)*math.pi*t/(2.0*self.H)) * math.sin((2.0*i-1.0)*math.pi*y/(2.0*self.H))
            u += tmp
            if (u == 0.0 or abs(tmp/u) < self.error_ratio or i >= max_iter_num):
                break
        u += -0.5 * self.density * self.bf * y*y / self.E + (self.p0 + self.density * self.bf * self.H) * y / self.E
        return u

    def stress(self, y, t):
        s = 0.0
        sign = 1.0
        i = 0
        max_iter_num = 100
        while (True):
            sign *= -1.0
            i += 1
            s_n = (2.0*i-1.0)*math.pi / (2.0*self.H) * 8.0*self.H * (2.0*math.pi*self.p0*i*sign - 2.0*self.density*self.bf*self.H - math.pi*self.p0*sign) / ((4.0*i*i - 4.0*i + 1.0) * (2.0*i - 1.0) * math.pi*math.pi*math.pi)
            tmp = s_n * math.cos(self.alpha*(2.0*i-1.0)*math.pi*t/(2.0*self.H)) * math.cos((2.0*i-1.0)*math.pi*y/(2.0*self.H))
            s += tmp
            if (s == 0.0 or abs(tmp/s) < self.error_ratio or i >= max_iter_num):
                break
        s += self.p0 + self.density * self.bf * (self.H - y)
        return s

if __name__ == "__main__":
    # parameters
    H = 1.0
    p0 = 1.0
    bf = 0.0
    E = 100.0
    density = 10.0
    t_len = 3.0
    data_num = 300
    # cal data
    bav = BarAxialVibration(H, p0, bf, E, density)
    t = np.zeros(data_num)
    u = np.zeros(data_num)
    s = np.zeros(data_num)
    t_inv = t_len / float(data_num)
    for i in range(data_num):
        t[i] = t_inv * i
        u[i] = bav.displacement(H, t[i])
        s[i] = bav.stress(H/2.0, t[i])
    # plot
    fig = plt.figure()
    plot1 = fig.subplots(1, 1)
    plot1.set_xlim(0.0, t_len)
    #plot1.plot(t, u)
    plot1.plot(t, s)
    plt.show()
