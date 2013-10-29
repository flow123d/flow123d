import math

def func_xyz(x,y,z):
    return ( x*y*z , )     # one value tuple

def func_circle(r,phi):
    return ( r * math.cos(phi), r * math.sin(phi) )
