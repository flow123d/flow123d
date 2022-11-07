import math

def func_xyz(x,y,z):
    return ( x*y*z , )     # one value tuple

def func_multi(x,y,z):
    return ( x*y*z )     # one value

def func_circle(r,phi):
    return ( r * math.cos(phi), r * math.sin(phi), 1 )
