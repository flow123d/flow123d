import numpy as np

# MANDEL'S 2D TEST CASE
# analytical solution according to:
# - R. Liu. Discontinuous Galerkin Finite Element Solution for Poromechanics, Ph.D. thesis (2004)
# - http://www.it.uu.se/edu/course/homepage/projektTDB/ht15/project08/Project08_Report.pdf
# (both contains some errors)
#
# other references for the description of the problem:
# - A. Verruijt. Theory and problems of poroelasticity. Delft, 2014
# - M. Ferronato et al. / Journal of Computational Physics 229 (2010) 4813–4830


# time unit
tunit = 1# 3600 # 1d=86400s

# half length of the domain in the horizontal direction
lengthA = 1

### elasticity parameters ###
# input lamb, G:
#lamb = 4e7
#G = 4e7
#K = lamb + 2/3*G
#E = 9*K*G / (3*K+G)
#poiss = (3*K - 2*G) / (2*(3*K+G))
#
#input E,poiss:
E = 1e0
poiss = 0.2
lamb = E*poiss/(1+poiss)/(1-2*poiss)
G = E / (2*(1+poiss))
K = E / (3*(1-2*poiss))
### end elasticity parameters ###

# hydraulic conductivity
hcond = 1e-0

# parameters for conversion of pressure to pressure head
# fluid density
rho = 1
# gravitational acceleration
gravity = 1

# external load on the top surface
force = 1e0

# Biot coefficient
alpha = 1

# storativity - can be calculated from compressibilities or set directly
# compressibilities (fluid, bulk)
#beta_f = 4.4e-4 #* 1e-6 # for lower storativity not converging in Flow123d
#beta_b = 0
# porosity
#por = 0.375
# storativity
#S = (por*beta_f + (alpha-por)*beta_b)
S = 1e-0





# Liu - nu_u (undrained Poisson's ratio)
poiss_u = (alpha**2*rho*gravity * (1+poiss) + 3*poiss*S*K) / (2*alpha**2*rho*gravity*(1+poiss) + 3*S*K)

# Liu - Skempton pore pressure coefficient (in Liu use the relation (4.1) for alpha, because (4.2) is wrong)
B = 3 * (poiss_u-poiss)/ ((1-2*poiss)*(1+poiss_u)*alpha)

# consolidation coefficient (Verruit) for Mandel's problem
cv = hcond / (rho * gravity) * (K+4/3*G) / (alpha**2*rho*gravity + S*(K+4/3*G))



# initial pressure, pressure head and displacement
p0 = 0.5 * force * alpha / (alpha**2 + S*(K+G/3))
h0 = p0 / (rho * gravity)



print("K: {:e}, G: {:e}".format(K,G))
print("E: {:e}, poiss: {}, poiss_u: {}".format(E,poiss,poiss_u))
print("lamb: {:e}, mu: {:e}".format(lamb,G))
print("force: {}".format(force))
print("k: {}".format(hcond))
print("S: {}".format(S))
print("B: {}".format(B))
print("cv: {}".format(cv))
print("p0: {}".format(p0))
print("h0: {}".format(h0))


def singular_points_func(idx, param):
    """Newton method for i-th solution of equation 'tan(x) - 2*param*x = 0'."""
    
    xk = np.pi * (idx - 1) + np.pi / 2.00001
    xkk = xk
    maxit = 1000
    tol = 1e-12
    f = np.tan(xk) - param * xk
    for i in range(maxit):
        c = np.cos(xk)
        xkk = xk - f * c * c / (1 - param * c * c)
        res1 = np.abs(xkk - xk) / np.abs(xkk)
        f = np.tan(xkk) - param * xkk
        res2 = np.abs(f)
        if res1 < tol and res2 < tol:
            break
        xk = xkk
    return xkk


singular_points = []
for i in range(1,200):
    singular_points.append(singular_points_func(i, (1-poiss)/(poiss_u-poiss)))  # according to Verruit param=2*eta, in Liu param=(1-poiss)/(poiss_u-poiss)


def analytic_sol(x, z, t):
    """ Analytical solution of pressure head and displacement."""

    s_p = 0
    s_u = 0
    s_v = 0
    for ksi in singular_points:
        cksi = np.cos(ksi)
        sksi = np.sin(ksi)
        a0 = sksi*(np.cos(ksi * x / lengthA) - cksi)
        a1 = ksi-sksi*cksi
        a2t = (ksi / lengthA) ** 2 * cv * t
        a2 = np.exp(-a2t)
        s_p = s_p + a0 / a1 * a2
        s_u = s_u - force*poiss_u/(2*G*lengthA)*sksi*cksi/a1*a2*x + force/G*cksi/a1*np.sin(ksi*x/lengthA)*a2
        s_v = s_v + sksi*cksi/a1*a2

    p = 2 * force * B * (1+poiss_u) / (3*lengthA) * s_p
    h = p / (rho * gravity)

    u = force*poiss/(2*G*lengthA)*x + s_u
    v = ( force*(1-poiss_u)/(G*lengthA)*s_v - force*(1-poiss)/(2*G*lengthA) )*z
    return h, u, v


def plot_presure_head():
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2)
    axs[0].set_ylabel('pressure [Pa]')
    axs[0].set_xlabel('x [m]')
    axs[1].set_ylabel('pressure [Pa]')
    axs[1].set_xlabel('t [s]')

    xax = np.arange(-lengthA, lengthA+0.05, 0.05)

    vals = np.zeros(xax.shape)
    t_vals = np.zeros(times.shape)
    ti = 0
    for t in times:
        for i in range(len(xax)):
            h, u, v = analytic_sol(xax[i], 0, t)
            vals[i] = h
            if abs(xax[i]) < 1e-6:
                t_vals[ti] = h
        axs[0].plot(xax/lengthA, vals, label="t={}".format(t))
        ti = ti + 1

    axs[0].legend()

    axs[1].plot(times, t_vals)
    print("pressure-over-time: {}".format(t_vals))

    plt.grid(color='lightgray', linestyle='--', linewidth=0.5)
    fig.suptitle("Analytical solution: pressure")
    #fig.tight_layout()  # otherwise the right y-label is slightly clipped
    
    plt.savefig("03_pressure_head.pdf")


def plot_displacement():
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2)
    axs[0].set_xlabel('x [m]')
    axs[0].set_ylabel('x_displacement [m]')
    axs[1].set_xlabel('x [m]')
    axs[1].set_ylabel('z_displacement [m]')

    xax = np.arange(-lengthA, lengthA+0.05, 0.05)

    u_vals = np.zeros(xax.shape)
    v_vals = np.zeros(xax.shape)
    for t in times:
        for i in range(len(xax)):
            h, u, v = analytic_sol(xax[i], lengthA, t)
            u_vals[i] = u
            v_vals[i] = v
        axs[0].plot(xax/lengthA, u_vals, label="t={}".format(t))
        axs[1].plot(xax/lengthA, v_vals, label="t={}".format(t))
       

    axs[0].legend()
    axs[1].legend()

    fig.suptitle("Analytical solution: displacement")
    #fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.savefig("03_displacement.pdf")


times = np.arange(0, 10, 1)
plot_presure_head()
plot_displacement()

### Paraview script

input0 = inputs[0]
ncells = input0.GetNumberOfCells()
time = input0.GetInformation().Get(vtk.vtkDataObject.DATA_TIME_STEP())

pressureArray = np.empty(ncells, dtype=np.float64)
displacementArray = np.empty(ncells, dtype=np.float64)
for i in range(ncells):
    cell = input0.GetCell(i)
    p1 = input0.GetPoint(cell.GetPointId(0))
    x1, y1, z1 = p1[:3]
    p2 = input0.GetPoint(cell.GetPointId(1))
    x2, y2, z2 = p2[:3]
    ttx = (x1+x2)/3
    tty = (y1+y2)/3
    ttz = (z1+z2)/3

    ttz = ttz + L
    h, u = analytic_sol(ttz, time * tunit)
    pressureArray[i] = h
    displacementArray[i] = u

output.CellData.append(pressureArray, "exact_pressure_head")
output.CellData.append(displacementArray, "exact_displacement")
