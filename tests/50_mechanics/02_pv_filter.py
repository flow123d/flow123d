import numpy as np
import math

# TERZAGHI'S 1D COLUMN TEST CASE
# DATA:
# according to
# https://www.researchgate.net/publication/222660034_A_fully_coupled_3-D_mixed_finite_element_model_of_Biot_consolidation/figures?lo=1
tunit = 1# 3600 # 1d=86400s

# elastic moduli
lamb = 4e7
G = 4e7
K = lamb + 2/3*G
E = 9*K*G / (3*K+G)
poiss = (3*K - 2*G) / (2*(3*K+G))

print("K: {:e}, G: {:e}".format(K,G))
print("E: {:e}, poiss: {}".format(E,poiss))
print("lamb: {:e}, mu: {:e}".format(lamb,G))

# hydraulic conductivity and density of fluid
k = 1e-5
rho = 1000

# porosity, Biot coefficient, compressibilities
por = 0.375
alpha = 1
beta_f = 4.4e-4 #* 1e-6 # for lower storativity not converging in Flow123d
beta_b = 0

# storage coefficient [Pa^{-1}]
S = (por*beta_f + (alpha-por)*beta_b)

# gravitational acceleration
g = 9.81

# column height
L = 0.25
# force on open end of column
q = 1e4

mv = 1.0/(K+4/3*G)
cv = k/(rho*g*(S+alpha*alpha*mv))
p0 = q * alpha*mv / (S+alpha*alpha*mv)
h0 = p0 / (rho*g)

u0 = p0*S/alpha

print("k: {}".format(k))
print("S: {} Pa^(-1) ~ {} m^(-1)".format(S, S*rho*g))
print("mv: {}".format(mv))
print("cv: {}".format(cv))
print("p0: {}".format(p0))
print("h0: {}".format(h0))
print("u0(top: z=0): {}".format(u0*(L)))

# def h_sum(z,t):
#     """ Analytical solution of pressure head."""
#     if t==0:
#         return h0
#
#     nk = 100
#     s = 0
#     for i in range(1,nk):
#         # Verruijt
#         a0 = (2*i-1) * np.pi / 2 / L
#         a1 = np.cos(a0*z)
#         a2 = np.exp(-a0**2 * cv*t)
#         s = s + (-1)**(i-1)/(2*i-1) * a1 * a2
#
#         # Ferronato 2010, range from 0, map z!
#         # a0 = (2*i+1) * np.pi/2 / L
#         # a1 = np.sin(a0*z)
#         # a2 = np.exp(-a0**2 * cv*t)
#         # s = s + 1/(2*i+1) * a1 * a2
#     s = p0*4/np.pi * s
#     s = s / (rho*g)
#     return s
#
#
# def u_sum(z, t):
#     """ Analytical solution of displacement."""
#     if t == 0:
#         return u0 * (L - z)
#
#     nk = 100
#     s = 0
#     for i in range(1, nk):
#         # Verruijt
#         a0 = (2 * i - 1) * np.pi / 2 / L
#         a1 = np.cos(a0 * z)
#         a2 = np.exp(-a0 ** 2 * cv * t)
#         s = s + 1/(2 * i - 1)**2 * a1 * a2
#
#         # Ferronato 2010, range from 0, map z!
#         # a0 = (2*i+1) * np.pi/2 / L
#         # a1 = np.cos(a0*z)
#         # a2 = np.exp(-a0**2 * cv*t)
#         # s = s + 1/(2*i+1)**2 * a1 * a2
#     # s = 1 - 8 / np.pi**2 * s
#     # u = L*p0 * (S/alpha + s*mv)
#
#     s = 8 / (np.pi ** 2) * s
#     u = mv*p0*((L-z) - L*s) + u0 * (L - z)
#     u = -u
#     return u

def analytic_sol_pressure(z, t):
    """ Analytical solution of pressure head and displacement."""
    if t == 0:
        h = h0
        return h

    if (cv*t/(L*L)) < 0.01:
        p = p0*math.erf( (L-z) / (2*math.sqrt(cv*t)) )
    else:
        nk = 1000
        s_h = 0

        for i in range(nk, 0, -1):
            # Verruijt
            a0 = (2 * i - 1) * np.pi / 2 / L
            a1 = np.cos(a0 * z)
            a2 = np.exp(-(a0 ** 2) * cv * t)
            s_h = s_h + ((-1) ** (i - 1)) / (2 * i - 1) * a1 * a2
    
            # Ferronato 2010, range from 0, map z!
            # a0 = (2*i+1) * np.pi/2 / L
            # a1 = np.sin(a0*z)
            # a2 = np.exp(-a0**2 * cv*t)
            # s = s + 1/(2*i+1) * a1 * a2
            #
            # a3 = np.cos(a0*z)
            # s = s + 1/(2*i+1)**2 * a3 * a2
        p = p0 * 4 / np.pi * s_h

    h = p / (rho * g)

    return h


def analytic_sol_displacement(z, t):
    """ Analytical solution of pressure head and displacement."""
    if t == 0:
        u = u0 * z
        return u

    s_u = 0
    nk = 1000

    for i in range(nk, 0, -1):
        # Verruijt
        a0 = (2 * i - 1) * np.pi / 2 / L
        a1_u = np.cos(a0 * (L-z))
        a2 = np.exp(-(a0 ** 2) * cv * t)
        s_u = s_u + 1 / (2 * i - 1) ** 2 * a1_u * a2
    
        # Ferronato 2010, range from 0, map z!
        # a0 = (2*i+1) * np.pi/2 / L
        # a1 = np.sin(a0*z)
        # a2 = np.exp(-a0**2 * cv*t)
        # s = s + 1/(2*i+1) * a1 * a2
        #
        # a3 = np.cos(a0*z)
        # s = s + 1/(2*i+1)**2 * a3 * a2

    u = 8 / (np.pi ** 2) * s_u
    u = mv * p0 * (z - L * u) + u0 * z
    return u


times = np.array([1., 2., 5., 10., 100.]) # time scale for higher storativity
zax = np.arange(0, L, L*0.01)

def plot_presure_head():
    import matplotlib.pyplot as plt
    fig, ax1 = plt.subplots()
    ax1.set_xlabel('pressure head [m]')
    ax1.set_ylabel('z [m]')

    vals = np.zeros(zax.shape)
    for t in times:
        for i in range(len(zax)):
            vals[i] = analytic_sol_pressure(zax[i], t)

        ax1.plot(vals, zax, label="t={}".format(t))

    # ax1.tick_params(axis='y')
    ax1.set_ylim([zax[0], zax[-1]] + (zax[-1]-zax[-2]))
    ax1.legend()

    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    # plt.show()
    plt.savefig("02_pressure_head.pdf")


def plot_displacement():
    import matplotlib.pyplot as plt
    fig, ax1 = plt.subplots()
    ax1.set_xlabel('displacement [m]')
    ax1.set_ylabel('z [m]')

    vals = np.zeros(zax.shape)
    for t in times:
        for i in range(len(zax)):
            vals[i] = analytic_sol_displacement(zax[i], t)
            # vals[i] = -u_sum(L-zax[i], t)

        ax1.plot(vals, zax, label="t={}".format(t))

    # ax1.tick_params(axis='y')
    ax1.set_ylim([zax[0], zax[-1]] + (zax[-1]-zax[-2]))
    ax1.legend()

    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    # plt.show()
    plt.savefig("02_displacement.pdf")


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
