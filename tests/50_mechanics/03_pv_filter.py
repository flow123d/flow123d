import numpy as np

# MANDEL'S 2D TEST CASE
# DATA:
# according to
# https://www.researchgate.net/publication/222660034_A_fully_coupled_3-D_mixed_finite_element_model_of_Biot_consolidation/figures?lo=1

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

#input E,poiss:
E = 1e1
poiss = 0.2
lamb = E*poiss/(1+poiss)/(1-2*poiss)
G = E / (2*(1+poiss))
K = E / (3*(1-2*poiss))
### end elasticity parameters ###

# hydraulic conductivity
hcond = 1e-1

# fluid density
rho = 1000

gravity = 9.81
force = 1e4

# Biot coefficient
alpha = 1
# compressibilities (fluid, bulk)
#beta_f = 4.4e-4 #* 1e-6 # for lower storativity not converging in Flow123d
#beta_b = 0
# porosity
#por = 0.375
# storativity
#S = (por*beta_f + (alpha-por)*beta_b)
S = 1e-1





# Verruit - eta
eta = (K + 4/3*G)/(2*G) * (alpha**2 + S*(K+G/3))/(alpha**2)

# Liu - nu_u (undrained Poisson's ratio)
#poiss_u = (1-poiss+2*eta*poiss)/(2*eta)
poiss_u = (alpha**2 * (1+poiss) + 3*poiss*S*K) / (2*alpha**2*(1+poiss) + 3*S*K)

# Liu - Skempton pore pressure coefficient (in Liu use the relation (4.1) for alpha, because (4.2) is wrong)
B = 3 * (poiss_u-poiss)/ ((1-2*poiss)*(1+poiss_u)*alpha)

# consolidation coefficient (Verruit) for Mandel's problem
cv = hcond / (rho * gravity) * (K+4/3*G) / (alpha**2 + S*(K+4/3*G))



# initial pressure, pressure head and displacement
p0 = 0.5 * force * alpha / (alpha**2 + S*(K+G/3))
h0 = p0 / (rho * gravity)
u0 = p0*S/alpha



print("K: {:e}, G: {:e}".format(K,G))
print("E: {:e}, poiss: {}, poiss_u: {}".format(E,poiss,poiss_u))
print("lamb: {:e}, mu: {:e}".format(lamb,G))
print("eta: {:e}".format(eta))
print("force: {}".format(force))
print("k: {}".format(hcond))
print("S: {}".format(S))
print("B: {}".format(B))
# print("mv: {}".format(mv))
print("cv: {}".format(cv))
print("p0: {}".format(p0))
print("h0: {}".format(h0))
#print("u0(top: z=0): {}".format(u0 * lengthA))


def singular_points_func(idx, param):
    """Newton method for i-th solution of equation 'tan(x) - 2*param*x = 0'."""
    
    xk = np.pi * (idx - 1) + np.pi / 2.00001 #2.001
    xkk = xk
    maxit = 1000
    tol = 1e-12
    f = np.tan(xk) - 2 * param * xk
    for i in range(maxit):
        c = np.cos(xk)
        xkk = xk - f / (1 / c / c - 2 * param)
        res1 = np.abs(xkk - xk) / np.abs(xkk)
        f = np.tan(xkk) - 2 * param * xkk
        res2 = np.abs(f)
        if res1 < tol and res2 < tol:
            break
        xk = xkk
#    print("idx: {}, xkk: {}".format(idx, xkk))
    return xkk


singular_points = []
for i in range(1,200):
    singular_points.append(singular_points_func(i, (1-poiss)/(poiss_u-poiss)/2))  # according to Verruit param=2*eta, in Liu param=(1-poiss)/(poiss_u-poiss)/2


def analytic_sol(x, z, t):
    """ Analytical solution of pressure head and displacement."""
    #if t == 0:
    #    h = h0
    #    u = u0 * z
    #    return h, u

    s_p = 0
    s_u = 0
    s_v = 0
    for ksi in singular_points:
        # Liu
        cksi = np.cos(ksi)
        sksi = np.sin(ksi)
        a0 = sksi*(np.cos(ksi * x / lengthA) - cksi)
        a1 = ksi-sksi*cksi
        a2t = (ksi / lengthA) ** 2 * cv * t
        a2 = np.exp(-a2t)
        s_p = s_p + a0 / a1 * a2
        s_u = s_u - force/(2*G*lengthA)*poiss_u*sksi*cksi/a1*a2*x + force/G*cksi/a1*np.sin(ksi*x/lengthA)*a2
        s_v = s_v + force*(1-poiss_u)/(G*lengthA)*sksi*cksi/a1*a2*z

    p = 2 * force * B * (1+poiss_u) / (3*lengthA) * s_p
    h = p / (rho * gravity)

    # u = 8 / (np.pi ** 2) * s_u
    # u = mv * p0 * (z - L * u) + u0 * z
    u = s_u + force*poiss/(2*G*lengthA)
    v = s_v + force*(1-poiss)/(2*G*lengthA)
    return h, u, v


def plot_presure_head():
    import matplotlib.pyplot as plt
    fig, ax1 = plt.subplots()
    ax1.set_ylabel('pressure [Pa]')
    ax1.set_xlabel('x [m]')

    times = np.arange(0, 1000, 100)*50
    # times = np.array([1, 10, 20, 60, 120, 600, 1200, 3600]) # time scale for lower storativity
    #times = np.array([0.05, 0.5, 2.5, 5, 10]) # time scale for higher storativity
    #times = np.array([1e3, 1e5, 8e5, 2e6, 3e6, 5e6, 1e7]) # time scale for higher storativity
    xax = np.arange(-lengthA, lengthA+0.05, 0.05)

    vals = np.zeros(xax.shape)
    for t in times:
        for i in range(len(xax)):
            h, u, v = analytic_sol(xax[i], 0, t)
            vals[i] = h

        ax1.plot(xax/lengthA, vals/h0, label="t={}".format(t))

    # ax1.tick_params(axis='y')
    # ax1.set_ylim([xax[0], xax[-1]] + (xax[-1]-xax[-2]))
    ax1.legend()

    plt.grid(color='lightgray', linestyle='--', linewidth=0.5)
    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    # plt.show()
    plt.savefig("03_pressure_head.pdf")


def plot_displacement():
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2)
    axs[0].set_xlabel('x [m]')
    axs[0].set_ylabel('x_displacement [m]')
    axs[1].set_xlabel('x [m]')
    axs[1].set_ylabel('z_displacement [m]')

    times = np.arange(0, 1000, 100)*50
    #times = np.array([1, 10, 20, 60, 120, 600, 1200, 3600]) # time scale for lower storativity
    #times = np.array([1e3, 1e5, 8e5, 2e6, 3e6, 5e6, 1e7]) # time scale for higher storativity
    #times = np.array([1e4, 1e5])  # time scale for higher storativity
    #zax = np.arange(0, 15, 0.1)
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
       

    # ax1.tick_params(axis='y')
    #ax1.set_ylim([xax[0], xax[-1]] + (xax[-1]-xax[-2]))
    axs[0].legend()
    axs[1].legend()

    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    # plt.show()
    #plt.yscale('log')
    plt.savefig("03_displacement.pdf")


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
