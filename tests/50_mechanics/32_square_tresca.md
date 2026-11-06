# Square with one separating Tresca fracture

Run [32_square_tresca.yaml](32_square_tresca.yaml). This case supplies numerical output for comparison with the derivation below.

## Geometry and material

The square is `(-0.5,0.5)^2`, split by the horizontal line `y=0`.
The mesh has 32 triangles and four fracture segments of length `0.25`.
All upper bulk triangles have smaller input IDs than all lower triangles, so the master side is consistently the upper side.
Fracture nodes run left to right.

Use unit bulk cross-section, `E=1`, `nu=0`, and hence `mu=0.5`.
The fracture has aperture `d=0.1`, minimum aperture `0.1`, `E_f=0.02`, `nu_f=0`, `mu_f=0.01`, and `fracture_sigma=1`. Contact therefore starts closed. Set compression `p=0.2` and Tresca strength `g=0.05`.
Disable gravity and poroelastic loading (`biot_alpha=0` and zero pressure).

## Analytical solution

Let `s=u_x(0+)-u_x(0-)`, let `D=u_x(0.5)-u_x(-0.5)`, and let `tau=sigma_xy`in both bulk halves.
The elastic fracture coupling contributes shear traction `k*s`, with `k=mu_f/d=0.1`.
The friction reaction has magnitude at most `g`.
For a frozen previous-step slip `s_old`, the equilibrium and friction laws are

```
D = s + tau/mu
tau = k*s + r
r in g*partial(abs(s-s_old))
```

For sticking, set `s=s_old` and check `abs(tau-k*s)<=g`.
For sliding in direction `eta=sign(s-s_old)`, use

```
s = (mu*D - eta*g)/(mu+k)
tau = mu*(D-s)
```

The displacement is affine on each bulk half:

```
u_x(x,y) = (tau/mu)*y + s/2       for y>0
u_x(x,y) = (tau/mu)*y - s/2       for y<0
u_y(x,y) = -p*y/E
u_z(x,y) = 0
u_fracture(x) = (0,0,0)
```

The bulk Cauchy stress is constant, with `sigma_xx=sigma_zz=0`, `sigma_yy=-p`, and `sigma_xy=sigma_yx=tau`.
These expressions are exactly representable by the P1 displacement space.
Displacement traces on the fracture must be read separately for the two bulk halves and the lower-dimensional fracture; their coordinates coincide.

## Loading sequence and expected output

Output is written at times 0, 1, and 2. The first reference slip is zero; the subsequent references are the preceding converged slips.

| Time | D | State | s | tau | contact_normal_stress | contact_tangential_stress |
| --- | --- | --- | --- | --- | --- | --- |
| 0 | 0.02 | Stick | 0 | 0.01 | -0.2 | (-0.01,0,0) |
| 1 | 0.2 | Positive sliding | 1/12 | 7/120 | -0.2 | (-0.05,0,0) |
| 2 | -0.2 | Reverse sliding | -1/12 | -7/120 | -0.2 | (0.05,0,0) |

Normal and tangential contact outputs above apply to the fracture cells; both fields are zero on bulk cells.
Tangential output is the reaction on the upper (master) body. It excludes the elastic fracture contribution `-k*s`.
For each length-0.25 segment the expected physical multipliers are `lambda=0.05` and `q=0.0025, 0.0125, -0.0125` respectively.
The aperture remains `0.1` in all three states.