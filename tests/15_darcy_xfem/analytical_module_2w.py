# script for ProgrammableFilter to compute
# analytical solution of square problem with fracture
#
# Script computes both - pressure values and velocity values



import math

class Well :
    # well
    sigma_w = 10.0
    P_w = 100.0
    rho_w = 0.03
    x_w = 3.33
    y_w = 3.33
    z_w = 0.0
    
    # conductivity of well
    k1 = 1.0
    # conductivity of aquifer
    k2 = 1.0
    
    # auxiliary variables to be filled
    a = 0.0
    
    def __init__(self, rho, sigma, P, x, y, k2, aa):
        self.rho_w = rho
        self.sigma_w = sigma
        self.P_w = P
        self.x_w = x
        self.y_w = y
        self.k2 = k2
        self.a = aa

    # distance from well w
    def r(self, x, y):
        xx = (x-self.x_w)
        yy = (y-self.y_w)
        return math.sqrt(xx*xx + yy*yy)
        
    ###############
    # fuction value for p2 - 2D pressure
    def p2(self, x, y):
        tr = self.r(x,y)
        if tr < self.rho_w:
            return self.P_w
        else:
            return self.a * math.log(tr)
    
    # fuction value for p2 - 2D pressure
    def vx2(self, x, y):
        tr = self.r(x,y)
        if tr < self.rho_w:
            return 0.0
        else:
            return - self.k2 * self.a * (x-self.x_w) / (tr*tr)
    
    # fuction value for p2 - 2D pressure
    def vy2(self, x, y):
        tr = self.r(x,y)
        if tr < self.rho_w:
            return 0.0
        else:
            return - self.k2 * self.a * (y-self.y_w) / (tr*tr)
    
class Parameters :
    
    k2 = 1e-3
    
    wells = []
     
    # source
    U = 0.0
    omg = 1
    
    def __init__(self):
        # wells
        self.wells.append(Well(0.03, 100, 150, 4.1, 4.3, self.k2, -52.2453721919971));
        self.wells.append(Well(0.03, 100, 100, 5.7, 5.9, self.k2, -40.6805832786121));
    
    
    def p2_fce_value(self, x, y):
        p_sin = self.U*math.sin(self.omg*x)
        p2 = 0
        for w in self.wells:
          p2 = p2 + w.p2(x,y)
          
        return p2 + p_sin
  
    def vx2_fce_value(self, x, y):
        vx2_sin = -self.k2*self.U * self.omg*math.cos(self.omg*x)
        vx2 = 0
        for w in self.wells:
          vx2 = vx2 + w.vx2(x,y)
          
        return vx2 + vx2_sin
  
    def vy2_fce_value(self, x, y):
        vy2 = 0
        for w in self.wells:
          vy2 = vy2 + w.vy2(x,y)
          
        return vy2


####################
# the module initialization (program is commented in this version, see the end of file)
      
params = Parameters()

########################
# provided functions

def velocity_2d(xx , yy, zz):
    vx2 = Parameters.vx2_fce_value(params,xx,yy)
    vy2 = Parameters.vy2_fce_value(params,xx,yy)
    return (vx2, vy2, 0.0)
    
def all_values_1d( xx , yy, zz):
    p1 = 0.0
    vx1 = 0.0
    return (p1, vx1, 0.0, 0.0, 1.0)
    
def all_values_2d( xx, yy, zz ):
    p2 = Parameters.p2_fce_value(params,xx,yy)
    vx2 = Parameters.vx2_fce_value(params,xx,yy)
    vy2 = Parameters.vy2_fce_value(params,xx,yy)
    return (p2, vx2, vy2, 0.0, 1.0)
    
def all_values_3d( xx, yy, zz ):
    p2 = Parameters.p2_fce_value(params,xx,yy)
    vx2 = Parameters.vx2_fce_value(params,xx,yy)
    vy2 = Parameters.vy2_fce_value(params,xx,yy)
    return (p2, vx2, vy2, 0.0, 1.0)
