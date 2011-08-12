# this test solves heat equation
#
#       \partial_t u - k \partial^2_x u = 0
#       
#       on (t,x) \in (0,\infty)x(0,1)
#
#       with u(0,x) = 0    and    u(t,0)=0 , u(t,1) = 100.0
#
#       k=0.02
#
#  solution in series form is:
#
#  u(t,x) = 100.0 x + 200.0/ \pi \sum_{n=1}^\infty  (-1)^n / n exp(-0.02 n^2 \pi^2 t ) sin(n \pi x)
#
#
#  folowing code is programmable filter for Paraview that computes the solution in barycenters of input 2D mesh
#

import math

class Parameters :

    #parametry
    k=0.02
    bc=100.0
    time=0.002

    alfa=[]
    un=[]
    B0 = 0
    u0 = 0

    N = 100000
    interval_size = 100
    tolerance = 1e-10
    

    # expresion for series for pressure on 2D
    def alfa_term(self,n):
      
      n_pi = math.pi * n
      if (n % 2) 
        # odd
        return -1/n * math.exp(-self.k*n_pi*n_pi*self.time )
      else   
        #even
        return 1/n * math.exp(-self.k*n_pi*n_pi*self.time )
      

 
    # all precalculations
    def precalculations(self) :
      #vypocet alfa, un
      for n in range (1,self.N+1,1):
        self.alfa.append(self.alfa_term(n))


    # fuction value for p2 - 2D pressure
    def p2_fce_value(self,x, y):
         series_sum = self.summarize(self.series_2d, x, y ) 
         return self.k * x + 2*self.bc / math.pi * series_sum
         

    # one term of series for 2D value
    def series_2d(self, n, tx, ty):    
        n_pi = (n+1)*math.pi
        return (self.alfa[n]*math.sin(n_pi*tx)
     

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

      
      
################################################
# the program


pdi = self.GetPolyDataInput()
pdo = self.GetOutput()
newData = vtk.vtkDoubleArray()
newData.SetName("analytical_pressure")
nc = pdi.GetNumberOfCells()

params = Parameters()
params.precalculations()

for i in range (nc):
      cell = pdi.GetCell(i)
      pb = cell.GetNumberOfPoints()
      if pb==3:
            p1 = pdi.GetPoint(cell.GetPointId(0))
            x1, y1, z1 = p1[:3]
            
            p2 = pdi.GetPoint(cell.GetPointId(1))
            x2, y2, z2 = p2[:3]
            
            p3 = pdi.GetPoint(cell.GetPointId(2))
            x3, y3, z3 = p3[:3]

            tx=1-math.fabs((x1+x2+x3)/3) * Parameters.x_scale
            ty=math.fabs((y1+y2+y3)/3) * Parameters.y_scale
            p2=Parameters.p2_fce_value(params, tx, ty )
            #p2=5
            newData.InsertNextValue(p2)
            pdo.GetCellData().AddArray(newData)





