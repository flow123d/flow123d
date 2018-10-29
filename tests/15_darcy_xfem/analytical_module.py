# script for ProgrammableFilter to compute
# analytical solution of square problem with fracture
#
# Script computes both - pressure values and velocity values



import math

class Parameters :

    # conductivity
    k1 = 1.0
    k2 = 1.0

    # wells
    sigma_w = 10.0
    P_w = 100.0
    #P_w = 40.0
    rho_w = 0.03    
    x_w = 3.33
    y_w = 3.33
    
    # Dirichlet BC
    #P2_d = 0.0
    P2_d = 40.0
    R = 10.0
    #P1 = 100.0

    # auxiliary variables to be filled
    a = 0.0
    b = 0.0
    
    # source
    U = 20.0
    omg = 1
    
    
    # all precalculations
    def precalculations(self):
        #self.a = (self.P_w - self.P2_d)/(math.log(self.rho_w/self.R))
        self.a = (self.P2_d - self.P_w)/(1.0/self.sigma_w/self.rho_w - math.log(self.rho_w/self.R))
        self.b = self.P2_d - self.a *math.log(self.R)

    # distance from well w
    def r(self, x, y):
        xx = (x-self.x_w)
        yy = (y-self.y_w)
        return math.sqrt(xx*xx + yy*yy)
        
    ###############  
    # fuction value for p2 - 2D pressure
    def p2_fce_value(self, x, y):
        tr = self.r(x,y)
        if tr < self.rho_w:
            return self.P_w
        else:
            return self.b + self.a * math.log(tr) 
    
    # fuction value for p2 - 2D pressure
    def vx2_fce_value(self, x, y):
        tr = self.r(x,y)  
        if tr < self.rho_w:
            return 0.0
        else:
            return - self.k2 * self.a * (x-self.x_w) / (tr*tr)
    
    # fuction value for p2 - 2D pressure
    def vy2_fce_value(self, x, y):
        tr = self.r(x,y)  
        if tr < self.rho_w:
            return 0.0
        else:
            return - self.k2 * self.a * (y-self.y_w) / (tr*tr)
       
    
    
    def p2_sin_source(self, x, y):
        p_sin = self.U*math.sin(self.omg*x)
        p2 = self.p2_fce_value(x,y) + p_sin
        return p2
  
    def vx2_sin_source(self, x, y):
        vx2_sin = -self.k2*self.U * self.omg*math.cos(self.omg*x)
        vx2 = self.vx2_fce_value(x,y) + vx2_sin
        return vx2
  
    def vy2_sin_source(self, x, y):
        return self.vy2_fce_value(x,y)
      
    #def p2_sin_source(self, x, y):
        #tr = self.r(x,y)
        #tt = 0.0;
        #if tr < self.rho_w:
            #tt = 0.0
        #else:
            #tt = (tr-self.rho_w) * (tr-self.rho_w)
        
        #p_sin = self.U*math.sin(self.omg*x) * tt
        #p2 = self.p2_fce_value(self,x,y) + p_sin
        #return p2
  
    #def vx2_sin_source(self, x, y):
        #tr = self.r(x,y)
        #tt = 0.0;
        #if tr < self.rho_w:
            #tt = (tr-self.rho_w)
        
        #vx2_sin = -self.k2*self.U*tt * (self.omg*math.cos(self.omg*x)*tt + 2*math.sin(self.omg*x)*(x-self.x)/tr)
        #vx2 = self.vx2_fce_value(self,x,y) + vx2_sin
        #return vx2
  
    #def vy2_sin_source(self, x, y):
        #return self.vy2_fce_value(self,x,y)
      
      
    # SIMPLE SOURCE to test convergence order: P: 2.0, U: 1.0
    
    #def p2_sin_source(self, x, y):
        #p2 = self.U*math.sin(self.omg*x)
        #return p2
  
    #def vx2_sin_source(self, x, y):
        #vx2 = -self.U*self.omg * math.cos(self.omg*x)
        #return vx2
  
    #def vy2_sin_source(self, x, y):
        #return 0.0
### end class Parameters    


####################
# the module initialization (program is commented in this version, see the end of file)
      
params = Parameters()
params.precalculations()

########################
# provided functions

def all_sin_source(xx , yy, zz):
    p2 = Parameters.p2_sin_source(params, xx, yy)
    vx2 = Parameters.vx2_sin_source(params, xx, yy)
    vy2 = Parameters.vy2_sin_source(params, xx, yy)
    return (p2, vx2, vy2, 0.0, 1.0)
  
#def all_sin_source(xx , yy, zz):
    #p2 = Parameters.p_sin_source(params, xx, yy)
    #vx2 = Parameters.vx_sin_source(params, xx, yy)
    #vy2 = Parameters.vy_sin_source(params, xx, yy)
    #return (p2, vx2, vy2, 0.0, 1.0)

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
    p3 = Parameters.p2_fce_value(params,xx,yy)
    vx3 = Parameters.vx2_fce_value(params,xx,yy)
    vy3 = Parameters.vy2_fce_value(params,xx,yy)
    vz3 = 0.0
    return (p3, vx3, vy3, vz3, 1.0) 
      
      
################################################
# the program ( Paraview ProgrammableFilter )


"""
pdi = self.GetPolyDataInput()
pdo = self.GetOutput()

nc = pdi.GetNumberOfCells()
PressureData = vtk.vtkDoubleArray()
PressureData.SetName("element_scalars")
PressureData.SetNumberOfValues(nc)

VelocityData = vtk.vtkDoubleArray()
VelocityData.SetName("element_vectors")
VelocityData.SetNumberOfComponents(3)
VelocityData.SetNumberOfTuples(nc)

params = Parameters()
params.precalculations()

p_idx=0
v_idx=0
for i in range (nc):
      cell = pdi.GetCell(i)
      pb = cell.GetNumberOfPoints()
      if pb==3: # triangle
            p1 = pdi.GetPoint(cell.GetPointId(0))
            x1, y1, z1 = p1[:3]
            
            p2 = pdi.GetPoint(cell.GetPointId(1))
            x2, y2, z2 = p2[:3]
            
            p3 = pdi.GetPoint(cell.GetPointId(2))
            x3, y3, z3 = p3[:3]
            
            ttx=(x1+x2+x3)/3
            tty=(y1+y2+y3)/3
            tx=1-math.fabs(ttx) * Parameters.x_scale
            ty=math.fabs(tty) * Parameters.y_scale
            
            p2=Parameters.p2_fce_value(params, tx, ty )
            vx2=Parameters.vx2_fce_value(params,tx,ty) * Parameters.x_scale * math.copysign(1,ttx)  *(-1)
            vy2=Parameters.vy2_fce_value(params,tx,ty) * Parameters.y_scale * math.copysign(1,tty)
            #p2=5
            PressureData.SetValue(p_idx, p2)
            p_idx+=1
            #VelocitData.InsertNextTupleValue([vx2,vy2,0.0])
            VelocityData.SetValue(v_idx, vx2)
            v_idx+=1
            VelocityData.SetValue(v_idx, vy2)
            v_idx+=1
            VelocityData.SetValue(v_idx, 0.0)
            v_idx+=1
            
      if pb==2: # line segment0
            p1 = pdi.GetPoint(cell.GetPointId(0))
            x1, y1, z1 = p1[:3]

            p2 = pdi.GetPoint(cell.GetPointId(1))
            x2, y2, z2 = p2[:3]

            ttx=(x1+x2)/2
            tx=1 - math.fabs(ttx) * Parameters.x_scale          
            p1 = params.p1_fce_value(tx , 0.0 )
            vx1 = 2*Parameters.v1_fce_value(params,tx,0.0) * Parameters.x_scale * math.copysign(1,ttx) *(-1)
            
            PressureData.SetValue(p_idx, p1)
            p_idx+=1
            VelocityData.SetValue(v_idx, vx1)
            v_idx+=1
            VelocityData.SetValue(v_idx, 0.0)
            v_idx+=1
            VelocityData.SetValue(v_idx, 0.0)
            v_idx+=1
            
pdo.GetCellData().AddArray(PressureData)
pdo.GetCellData().AddArray(VelocityData)
"""