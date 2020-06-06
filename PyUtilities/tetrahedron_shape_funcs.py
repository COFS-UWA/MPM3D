import numpy as np

class Point:
    def __init__(self):
        self.x = 0.0
        self.y = 0.0
        self.z = 0.0
    
    def __init__(self, x, y ,z):
        self.x = x
        self.y = y
        self.z = z

    def div_num(self, num):
        self.x /= num
        self.y /= num
        self.z /= num
        return self
    
    def __str__(self):
        return "(%f, %f, %f)\n" % (self.x, self.y, self.z)

# form vector p2p1
def vec(p2, p1):
    return Point(p2.x-p1.x, p2.y-p1.y, p2.z-p1.z)

def dot_prod(v1, v2):
    return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z

# v1 x v2
def cross_prod(v1, v2):
    x = v1.y * v2.z - v2.y * v1.z
    y = v1.z * v2.x - v2.z * v1.x
    z = v1.x * v2.y - v2.x * v1.y
    return Point(x, y, z)

def cal_vol(p1, p2, p3, p4):
    v21 = vec(p2, p1)
    v31 = vec(p3, p1)
    v41 = vec(p4, p1)
    return dot_prod(cross_prod(v21, v31), v41) / 6.0

if __name__ == "__main__":
    p  = Point(0.3, 0.3, 0.3)
    p1 = Point(0.0, 0.0, 0.0)
    p2 = Point(1.0, 0.0, 0.0)
    p3 = Point(0.0, 1.0, 0.0)
    p4 = Point(0.0, 0.0, 1.0)
    vol = cal_vol(p1, p2, p3, p4)
    print("vol: %f\n" % vol)
    # N1, N2, N3, N4
    N1 = cal_vol(p2, p4, p3, p) / vol
    N2 = cal_vol(p1, p3, p4, p) / vol
    N3 = cal_vol(p1, p4, p2, p) / vol
    N4 = cal_vol(p1, p2, p3, p) / vol
    print("N1: %f\nN2: %f\nN3: %f\nN4: %f\n" % (N1, N2, N3, N4))
    # derivatives
    v21 = vec(p2, p1)
    v31 = vec(p3, p1)
    v41 = vec(p4, p1)
    v42 = vec(p4, p2)
    v32 = vec(p3, p2)
    # N1 derivatives
    dN1 = cross_prod(v42, v32).div_num(6.0*vol)
    print("dN1_dx: %f\ndN1_dy: %f\ndN1_dz: %f\n" % (dN1.x, dN1.y, dN1.z))
    # N2 derivatives
    dN2 = cross_prod(v31, v41).div_num(6.0*vol)
    print("dN2_dx: %f\ndN2_dy: %f\ndN2_dz: %f\n" % (dN2.x, dN2.y, dN2.z))
    # N3 derivatives
    dN3 = cross_prod(v41, v21).div_num(6.0*vol)
    print("dN3_dx: %f\ndN3_dy: %f\ndN3_dz: %f\n" % (dN3.x, dN3.y, dN3.z))
    # N4 derivatives
    dN4 = cross_prod(v21, v31).div_num(6.0*vol)
    print("dN4_dx: %f\ndN4_dy: %f\ndN4_dz: %f\n" % (dN4.x, dN4.y, dN4.z))
