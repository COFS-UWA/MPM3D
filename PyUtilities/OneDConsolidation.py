import math
import numpy as np
import matplotlib.pyplot as plt

class OneDConsolidation:
    """
    z = 0, free flow boundary condition
    z = H, impermeable boundary condition
    Parameters:
        1. Cv, coefficient of consolidation;
        2. Es, one dimensional compressive modulus
        3. u0, initial pore pressure;
        4. H, depth of soil;
        5. error_ratio, used to control the calculation precision.
    """
    def __init__(self, Cv, Es, u0, H, error_ratio = 1.0e-3):
        self.Cv = Cv
        self.Es = Es
        self.u0 = u0
        self.H = H
        # Final settlement
        self.dH_final = -H * u0 / Es
        self.error_ratio = error_ratio
    
    def calPorePressure(self, t, z):
        Tv = self.Cv * t / (self.H * self.H)
        p = 0.0
        z = z / self.H
        i = 0
        while True:
            M = (2*i+1) * math.pi / 2.0
            inc = 2.0/M * math.sin(M*z) * math.exp(-M*M*Tv)
            p += inc
            i += 1
            if abs(inc) < self.error_ratio:
                break
        if (p > 1.0): p = 1.0
        p *= self.u0
        return p
    
    def calSettlement(self, t):
        Tv = self.Cv * t / (self.H * self.H)
        dH = 0.0
        i = 0
        while True:
            M = (2*i+1) * math.pi / 2.0
            inc = 2.0/(M*M) * math.exp(-M*M*Tv)
            dH += inc
            i += 1
            if abs(inc) < self.error_ratio:
                break
        dH = self.dH_final * (1.0 - dH)
        return dH
    
if __name__ == "__main__":
    Es = 40.0e6
    kv = 1.0e-5
    miu = 1.0 # dynamic viscosity
    Cv = kv * Es / miu
    u0 = 40.0e3
    H = 10.0
    con_res = OneDConsolidation(Cv, Es, u0, H)
    
    fig = plt.figure()
    plot1 = fig.subplots(1, 1)
    plot1.set_title('Settlement - Time relation')
    plot1.set_xlabel('Time')
    plot1.set_ylabel('Settlement')
    
    data_num = 100
    t_list = np.zeros(data_num)
    p_list = np.zeros(data_num)
    u_list = np.zeros(data_num)
    for i in range(data_num):
        t_list[i] = 0.01 * float(i)
        p_list[i] = con_res.calPorePressure(t_list[i], 10.0)
        u_list[i] = con_res.calSettlement(t_list[i])
    
    plot1.set_xlim([t_list[0], t_list[data_num-1]])
    
    plot1.plot(t_list, p_list, 'k--')
    #plot1.plot(t_list, u_list, 'k--')
    
    plt.show()
