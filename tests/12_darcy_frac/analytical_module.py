# script for ProgrammableFilter to compute
# analytical solution of square problem with fracture
#
# Script computes both - pressure values and velocity values



import math

class Parameters :

    # parametry
    # only half of the fracture means half of the conductivity
    k1=1.0 *0.5
    k2=5.0

    # sigma has to be half of the sigma on 
    # the full square <-1,1> x <-1,1>
    sigma=100.0 *0.5
    P2=10.0
    P1=5.0

    x_scale = 1.0
    y_scale = 1.0

    alfa=[]
    un=[]
    B0 = 0
    u0 = 0

    N = 100000
    interval_size = 20
    tolerance = 1e-10
    

    # expresion for series for pressure on 2D
    def alfa_term(self,n):
      
      n_pi = math.pi * n
      k1_n_pi_n_pi = self.k1 * n_pi * n_pi
      
      a1=self.sigma*self.k2
      a2=self.k2 * n_pi * (self.sigma+ k1_n_pi_n_pi)
      a3=( 1 + math.exp(-2 * n_pi) )/(2)
      a4=self.sigma * k1_n_pi_n_pi
      a5=( 1 - math.exp(-2 * n_pi) )/(2)
      return ((a1/((a2*a3)+(a4*a5))))

    # expression for series for pressure on 1D    
    def u_term(self, n):
      
      u1=self.alfa[n-1] * self.sigma * ( 1 - math.exp(-2 * math.pi * n) )
      u2=(2*self.sigma)+(2*self.k1*math.pi*math.pi*n*n)
      return (u1/u2)

    # expression of absolute term B0
    def B0_term(self):
    
      b1=(self.P2-self.P1)*self.sigma*math.sqrt(self.k1)
      b2=math.sqrt(self.sigma) * self.k2 * \
           (math.cosh(math.sqrt(self.sigma/self.k1)) \
           /math.sinh(math.sqrt(self.sigma/self.k1)))
      b3=self.sigma*math.sqrt(self.k1)
      b4=2*self.sigma*math.sqrt(self.k1)

      b5=0
      for i in range (self.N-1,-1,-1):
        b5=b5+self.un[i]

      return (b1)/( b2+b3+(b4*b5) )
 
    # all precalculations
    def precalculations(self) :
      #vypocet alfa, un
      for n in range (1,self.N+1,1):
        self.alfa.append(self.alfa_term(n))
        self.un.append(self.u_term(n))

      self.B0=self.B0_term()
      
      #vypocet u0
      u0t=math.sqrt(self.k1 * self.sigma)*math.sinh(math.sqrt(self.sigma/self.k1))
      self.u0=(self.k2*self.B0)/(u0t)

    ###############  
    # fuction value for p2 - 2D pressure
    def p2_fce_value(self,x, y):
         series_sum = self.summarize(self.series_2d, x, y ) 
         return self.P2 - self.B0 + (self.B0 * y) - (2 * self.B0 * series_sum)

    # one term of series for 2D pressure value
    def series_2d(self, n, tx, ty):    
        n_pi = (n+1)*math.pi
        tmp1=math.cos(n_pi * tx)
        tmp2=math.exp(-n_pi * ty) * ( 1 - math.exp( -2 * n_pi * (1-ty) ) )
        return (self.alfa[n]*tmp1*tmp2)/2.0
        
    ############# 
    # fuction value for p1 - 1D pressure
    def p1_fce_value(self, x, y):
         series_sum = self.summarize(self.series_1d, x, y ) 
         pa=self.u0*math.cosh((1-x)*math.sqrt(self.sigma/self.k1))
         return self.P2-self.B0-pa-2*self.B0*series_sum
    # one term of series for 1D pressure value 
    def series_1d(self, n, tx, ty):     
        return (self.un[n]*math.cos(math.pi*(n+1)*tx))
        
    # fuction value for vx2 - x component of 2D velocity 
    def vx2_fce_value(self,x, y):
         series_sum = self.summarize(self.series_vx_2d, x, y ) 
         dp_dx= - (2 * self.B0 * series_sum)
         return -self.k2*dp_dx
    # one term of series for 2D vx value
    def series_vx_2d(self, n, tx, ty):    
        n_pi = (n+1)*math.pi
        tmp1=-n_pi*math.sin(n_pi * tx)
        tmp2=math.exp(-n_pi * ty) * ( 1 - math.exp( -2 * n_pi * (1-ty) ) )
        return (self.alfa[n]*tmp1*tmp2)/2.0

    # fuction value for vy2 - y component of 2D velocity 
    def vy2_fce_value(self,x, y):
         series_sum = self.summarize(self.series_vy_2d, x, y ) 
         dp_dy=self.B0  - (2 * self.B0 * series_sum)
         return -self.k2*dp_dy
    # one term of series for 2D value
    def series_vy_2d(self, n, tx, ty):    
        n_pi = (n+1)*math.pi
        tmp1=math.cos(n_pi * tx)
        tmp2=-n_pi*math.exp(-n_pi * ty) * ( 1 + math.exp( -2 * n_pi * (1-ty) ) )
        return (self.alfa[n]*tmp1*tmp2)/2.0

    # fuction value for v1 - 1D velocity
    def v1_fce_value(self, x, y):
         series_sum = self.summarize(self.series_v1_1d, x, y ) 
         pa=-self.u0*(math.sqrt(self.sigma/self.k1))*math.sinh((1-x)*math.sqrt(self.sigma/self.k1))
         dp_dx=-pa-2*self.B0*series_sum
         return -self.k1*dp_dx
    # one term of series for 1D value 
    def series_v1_1d(self, n, tx, ty):     
        return -(self.un[n]*math.pi*(n+1)*math.sin(math.pi*(n+1)*tx))
    
    # sum the series with terms given by function "series"
    def summarize( self, series, x, y):
        local_sum=[]
            
        interval_begin=0        
        for interval_end in range (self.interval_size, self.N, self.interval_size):
            # sum over the interval backward  
            interval_sum = 0.0
            for n in range (interval_end-1,interval_begin-1,-1):
                interval_sum += series( n, x, y)
                    
            if (math.fabs(interval_sum) < self.tolerance) :
                break    
            local_sum.append( interval_sum )
            interval_begin +=self.interval_size
     
        total_sum = 0.0
        for number in reversed(local_sum) :
            total_sum += number
        return total_sum    
            
### end class Parameters    


####################
# the module initialization (program is commented in this version, see the end of file)
      
params = Parameters()
params.precalculations()

########################
# provided functions

def all_values_1d( xx , yy, zz):
    tx=1-math.fabs( xx) * Parameters.x_scale
    p1 = params.p1_fce_value(tx , 0.0 )
    vx1 = 2*Parameters.v1_fce_value(params,tx,0.0) * Parameters.x_scale * math.copysign(1,xx) *(-1)
    return (p1, vx1, 0.0, 0.0, 1.0)
        
def all_values_2d( xx, yy, zz ):
    tx=1-math.fabs( xx ) * Parameters.x_scale
    ty=math.fabs( yy ) * Parameters.y_scale
    p2=Parameters.p2_fce_value(params, tx, ty )
    vx2=Parameters.vx2_fce_value(params,tx,ty) * Parameters.x_scale * math.copysign(1,xx)  *(-1)
    vy2=Parameters.vy2_fce_value(params,tx,ty) * Parameters.y_scale * math.copysign(1,yy)
    return (p2, vx2, vy2, 0.0, 1.0) 
    

      
      
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