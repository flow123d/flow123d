import numpy as np

# MANDEL'S 2D TEST CASE
# DATA:
# according to
# https://www.researchgate.net/publication/222660034_A_fully_coupled_3-D_mixed_finite_element_model_of_Biot_consolidation/figures?lo=1
tunit = 1# 3600 # 1d=86400s

lamb = 4e7
G = 4e7
K = lamb + 2/3*G
E = 9*K*G / (3*K+G)
poiss = (3*K - 2*G) / (2*(3*K+G))

# mv = 1.0/(K+4/3*G)
nu = (K + 4/3*G)/(2*G)

print("K: {:e}, G: {:e}".format(K,G))
print("E: {:e}, poiss: {}".format(E,poiss))
print("lamb: {:e}, mu: {:e}".format(lamb,G))
print("nu: {:e}".format(nu))

hcond = 1e-5
rho = 1000

por = 0.375
alpha = 1
beta_f = 4.4e-4 #* 1e-6 # for lower storativity not converging in Flow123d
beta_b = 0
S = (por*beta_f + (alpha-por)*beta_b)


gravity = 9.81
force = 1e4

cv = hcond / (rho * gravity) * (K+4/3*G)

p0 = 0.5 * force
h0 = p0 / (rho * gravity)

u0 = p0*S/alpha

lengthA = 1

print("force: {}".format(force))
print("k: {}".format(hcond))
print("S: {}".format(S))
# print("mv: {}".format(mv))
print("cv: {}".format(cv))
print("p0: {}".format(p0))
print("h0: {}".format(h0))
print("u0(top: z=0): {}".format(u0 * lengthA))


def singular_points_func(idx, param):
    """Newton method for i-th solution of equation 'tan(x) - 2*nu*x'."""
    xk = np.pi * (idx - 1) + np.pi / 2.001
    xkk = xk
    maxit = 1000
    tol = 1e-12
    f = np.tan(xk) - 2 * param * xk
    for idx in range(maxit):
        c = np.cos(xk)
        xkk = xk - f / (1 / c / c - 2 * param)
        res1 = np.abs(xkk - xk) / np.abs(xkk)
        f = np.tan(xkk) - 2 * param * xkk
        res2 = np.abs(f)
        if res1 < tol and res2 < tol:
            break
        xk = xkk
    return xkk


singular_points = []
for i in range(1, 100):
    singular_points.append(singular_points_func(i, nu))


def analytic_sol(x, z, t):
    """ Analytical solution of pressure head and displacement."""
    if t == 0:
        h = h0
        u = u0 * z
        return h, u

    nk = 100
    s_h = 0
    # s_u = 0
    for ksi in singular_points:
        # Verruijt
        cksi = np.cos(ksi)
        a0 = cksi*(np.cos(ksi * x / lengthA) - cksi)
        a1 = 1-2*nu*cksi**2
        a2t = (ksi / lengthA) ** 2 * cv * t
        a2 = np.exp(-a2t)
        s_h = s_h + a0 / a1 * a2
        # s_u = s_u + 1 / (2 * i - 1) ** 2 * a1_u * a2

    p = p0 * 2 * nu * s_h * 2
    h = p / (rho * gravity)

    # u = 8 / (np.pi ** 2) * s_u
    # u = mv * p0 * (z - L * u) + u0 * z
    u = u0
    return h, u


def plot_presure_head():
    import matplotlib.pyplot as plt
    fig, ax1 = plt.subplots()
    ax1.set_ylabel('pressure [Pa]')
    ax1.set_xlabel('x [m]')

    # times = np.arange(60, 3600, 60)
    # times = np.array([1, 10, 20, 60, 120, 600, 1200, 3600]) # time scale for lower storativity
    times = np.array([0.05, 0.5, 2.5, 5, 10]) # time scale for higher storativity
    xax = np.arange(-lengthA, lengthA+0.05, 0.05)

    vals = np.zeros(xax.shape)
    for t in times:
        for i in range(len(xax)):
            h, u = analytic_sol(xax[i], 0, t)
            vals[i] = h

        ax1.plot(xax/lengthA, vals/h0, label="t={}".format(t))

    # ax1.tick_params(axis='y')
    # ax1.set_ylim([xax[0], xax[-1]] + (xax[-1]-xax[-2]))
    ax1.legend()

    plt.grid(color='lightgray', linestyle='--', linewidth=0.5)
    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    # plt.show()
    plt.savefig("05_pressure_head.pdf")


def plot_displacement():
    import matplotlib.pyplot as plt
    fig, ax1 = plt.subplots()
    ax1.set_xlabel('displacement [m]')
    ax1.set_ylabel('z [m]')

    # times = np.arange(60, 3600, 60)
    # times = np.array([1, 10, 20, 60, 120, 600, 1200, 3600]) # time scale for lower storativity
    # times = np.array([1e4, 1e5, 1e6, 1e7, 5e7]) # time scale for higher storativity
    times = np.array([1e4, 1e5])  # time scale for higher storativity
    zax = np.arange(0, 15, 0.1)

    vals = np.zeros(zax.shape)
    for t in times:
        for i in range(len(zax)):
            h, u = analytic_sol(zax[i], t)
            vals[i] = u
            # vals[i] = -u_sum(L-zax[i], t)

        ax1.plot(vals, zax, label="t={}".format(t))

    # ax1.tick_params(axis='y')
    ax1.set_ylim([zax[0], zax[-1]] + (zax[-1]-zax[-2]))
    ax1.legend()

    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    # plt.show()
    plt.savefig("05_displacement.pdf")


plot_presure_head()
# plot_displacement()

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
