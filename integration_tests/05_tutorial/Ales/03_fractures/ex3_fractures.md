
# FLOW123D EXAMPLES
## 1. Example 3 – fractures
In Flow123D domain interaction with fractures is possible to implemented. This example was created for evaluation of inluence of individual processes (difusion, linear sorption, dual-porosity) between domain interaction in transport. The background of this tast is movement of contaminant mass from deep repository along fractures. The output mass from fracture is evaluated.
 
### 1.1    Geometry and boundary conditions
The simulation area is in the range 100 × 100 m with 1 m thicness. The area is cut by two flow fractures and from them two blind fractures are separeted. The cross-section of fractures is 0.01 m. Geometry and mesh is depicted in Figure 1.

![alt text](http://bacula.nti.tul.cz/~ales.balvin/ex3_fractures/geometry)

Figure 1: Geometry of simulation area.

![alt text](http://bacula.nti.tul.cz/~ales.balvin/ex3_fractures/mesh)

Figure 2: Mesh of simulation area. 

### Mesh regions - Union in code
The "!Union" instruction in the mesh section was used for simplification and shortening yaml file. By this function is possible to adapt mesh boundary condition ( Control file: ex3_fractures_diffusion.yaml):

  mesh:
    mesh_file: ./input/ex3_mesh.msh
    regions:
      - !Union
        name: flow_fractures
        regions:
          - flow_fracture1
          - flow_fracture2
      - !Union
        name: blind_fractures
        regions:
          - blind_fracture1
          - blind_fracture2
      - !Union
        name: BC_right
        regions:
          - .right
          - .right_points
      - !Union
        name: BC_left
        regions:
          - .left
          - .left_points
    

### 1.1.2 Hydraulic boundary condition

Four variants of the model was considered which differed by hydraulic conductivity and present of the blind fractures (Table 1.) in ex3_fractures_diffusion.yaml: 

    input_fields:
    - region: rock
      conductivity: 0.0315
      cross_section: 1
    - region: flow_fractures
      conductivity: 31.5
      cross_section: 0.01
    - region: blind_fractures
      conductivity: 3.15
      cross_section: 0.01


Table 1.
Vaiants    |Blind fractures | K rock     [m/s] | K fracture [m/s] | n rock     | n fracture
-----------|----------------|------------------|------------------|------------|--------
Variant 1 |     Yes        |1e-9 (0.0315 m/yr)|1e-6 (31.5 m/yr)  |    0.005   |   0.1
Variant 2 |      Yes       |     1e-11        |    1e-6          |    0.005   |   0.1
Variant 3 |      No        |     1e-9         |    1e-6          |    0.005   |   0.1
Variant 4 |      No        |     1e-11        |    1e-6          |    0.005   |   0.1

When the blind fractures are present in model, that their parameters are identical with rock. On the oposite variants, their parameters are in accordance with flow-fractures.
    
Two dirichlet boundary conditions were definated for flux - piezometric head 0.1 m on the left side and 0 m on the right site (ex3_fractures_diffusion.yaml): 

    - region: BC_left
      bc_type: dirichlet
      bc_piezo_head: 0.1
    - region: BC_right
      bc_type: dirichlet
      bc_piezo_head: 0

Other sides were nonpermeable. Boundary condition was considered with this values because of the filtration flux of flux-fractures was in order of magnitude 1e-9 m/s (it is 0.1 m/yr).

### 1.1.3 Transport boundary condition

The whole advection-diffusion solution is used for solving this task problem:

    solute_equation: !Coupling_OperatorSplitting
      transport: !Solute_AdvectionDiffusion_DG

The transport boundary condition was assumed on the fracture in Gauss curve, when the middle value was prescribed 2000 years and standart deviation 700 years. 
Equation for Gauss curve is in shape:

![alt text](http://bacula.nti.tul.cz/~ales.balvin/ex3_fractures/eq)

For diffusion, the Gauss curve is prescribed in the yaml code as (ex3_fractures_flux_diffusion.yaml):

      - region: .left_0
        bc_type: dirichlet
        bc_conc: !FieldFormula
          value: (2.84959e-5)*exp(-(t-2000)^2/980000)

That means that during the simulation time is relased 0.05 kg/m3 flow water. Maximum concentration of realesed water is 0.028 g/m3. Mean value correspond with rael value of release of isotopes of deep repository.

The parameters of mechanical dispersion are set on the some value for all variants of model (longitudinaly dispersivity = 5 m and transversaly dispersivity 0.5 m). The diffusion is evaluated according to variants (ex3_fractures_diffusion.yaml): 

      input_fields:
      - region: rock
        init_conc: 0
        porosity: 0.005
        diff_m: 0.0369
        disp_l: 5
        disp_t: 0.5
      - region: flow_fractures
        init_conc: 0
        porosity: 0.1
        diff_m: 0.0369
        disp_l: 5
        disp_t: 0.5
      - region: blind_fractures
        init_conc: 0
        porosity: 0.1
        diff_m: 0.0369
        disp_l: 5
        disp_t: 0.5
    
## 1.2    Diffusion

In this part the coefficient of molekular diffusion was changed in range from 1e-14 to 1e-11 m2/s. In Flow123d the value of coefficient is prescribed in control file after division of multiple of porosity by cude root of porosity. This value was the some at rock and fractures, that means that in reality the value on the fractures is 1.5 higher then on the rock (summary in Table 2.)
     
Tab. 2: Coefficient of molekular diffusion prescribed in Flow123d

Coefficient of molekular diffusion Flow123D [m2/s] | Coefficient of efective molekular diffusion – rock [m2/s] |  Coefficient of efective molekular diffusion – fracture [m2/s] 
-------------|-----------|-----------
1.16961e-8   |   1e-11   | 5.43e-10  
1.16961e-9   |   1e-12   | 5.43e-11
1.16961e-10  |   1e-13   | 5.43e-12
1.16961e-11  |   1e-14   | 5.43e-13

### 1.2.1  Diffusion - Results

Time dependent of mass balance from fracture is depicted in Figure 3 and maximum of mass is in Tables 3-6.
In all model variants maximum of mass flux decreases with increases of molekular diffusion coeficient. Presence of blind fractures has minimal influence to results (decrease of maximal flux with higher diffusion was a little higher then in variants without blind fractures). Overall the influence of diffusion is more higher in variants with lower conductivity. While in variant with minimal difusion (1e-14) is maximal mass flux from fractures in variants with low conductivity more then five times higher then variants with higer conductivity. In variant with maximal difusion is one and half. 

![alt text](http://bacula.nti.tul.cz/~ales.balvin/ex3_fractures/mass_flux_fracture.png)

Figure 3:  Maximum of mass flux from fracture, (a) higher conductivity for rock and with blind fractures
(b) higher conductivity for rock and with blind fractures, (c) lower conductivity for rock and with blind fractures,
(d) higher conductivity for rock and without blind fractures.

Diffusion coeficient [m2/s] | Maximum [kg] | Change again variant with Diff=1e-14 [%]
-------------|-----------|-----------
1e-14        | 1.08e-9   |    0
1e-13        | 1.03e-9   | -4.6
1e-12        | 8,19e-10  |-24,2
1e-11        | 7,79e-10  |-27,9
Tab. 3: Maximum of mass flux from fracture - higher conductivity for rock and with blind fractures

Diffusion coeficient [m2/s] | Maximum [kg] | Change again variant with Diff=1e-14 [%]
-------------|-----------|-----------
1e-14        | 5,81e-9   |0
1e-13        | 3,14e-9   |-46
1e-12        | 5,27e-10  |-90,9
1e-11        | 4,63e-10  |-92
Tab. 4: Maximum of mass flux from fracture - lower conductivity for rock and with blind fractures

Diffusion coeficient [m2/s] | Maximum [kg] | Change again variant with Diff=1e-14 [%]
-------------|-----------|-----------
1e-14        | 1,08e-9        | 0
1e-13        | 1,02e-9        | -5,6
1e-12        | 7,71e-10        | -28,6
1e-11        | 7,35e-10        | -31,9
Tab. 5: Maximum of mass flux from fracture - higher conductivity for rock and without blind fractures

Diffusion coeficient [m2/s] | Maximum [kg] | Change again variant with Diff=1e-14 [%]
-------------|-----------|-----------
1e-14        |7,05e-9        |0
1e-13        |3,57e-9        |-49,4
1e-12        |5,14e-10        |-92,7
1e-11        |4,47e-10        |-93,7
Tab. 6: Maximum of mass flux from fracture - higher conductivity for rock and without blind fractures

## 1.3    Sorption
Flow123d allow computing sorption with linear, Langmuir and Freundlich isoterm. In this part, five elements was computing with various isotherm. See in yaml file ("ex3_sorption.yaml"):

    reaction_term: !Sorption
      substances:
        - A
        - As-lin-limit
        - As-lin
        - As-lang
        - As-freund
      solvent_density: 1.0
      solubility:
        - 1.0
        - 1.0
        - 0.1
        - 1.0
        - 1.0
      input_fields:
        - region: all_rock
          init_conc_solid:
            - 0
            - 0
            - 0.5
            - 0
            - 0
          rock_density: 1.0
          sorption_type:
            - linear             # reach solubility limit
            - linear
            - none            # no sorption - trying to switch off
            - langmuir
            - freundlich
          isotherm_mult: 0.6
          isotherm_other: 0.4


![alt text](http://bacula.nti.tul.cz/~ales.balvin/ex3_fractures/ex3_resuts_sorption_curves)
Figure 4: 


## 1.4    Dual-porosity
Dual porosity substitutes blind fractures in this task. The dual-porosity parameters ("diffusion_rate_immobile": 5.64742e-06) was calibrated for identical results with the model with the blind fractures. Other settings of transport is identical with previous diffusion model (without diffusion and dispersion for "transport: !Solute_Advection_FV" variant).
Dual porosity was compute with numerical difusion "transport: !Solute_Advection_FV" (for "!Solute_Advection_FV" is necessery calibration other parameters). In this section a new section "reaction_term: !DualPorosity" for dual porosity is set (ex3_dual_porosity.yaml):

    reaction_term: !DualPorosity
      input_fields:
        - region: rock
          init_conc_immobile:
            - 0
        - region: flow_fractures
          diffusion_rate_immobile: 5.64742e-06
          porosity_immobile: 0.01
          init_conc_immobile:
            - 0
        - region: blind_fractures
          init_conc_immobile:
            - 0

### 1.4.1    Dual-porosity calibration and comparison with difusion
Results of calibration of model with dual porosity and model with flow in blin_fractures is depicted in Fig. 4.  


![alt text](http://bacula.nti.tul.cz/~ales.balvin/ex3_fractures/ex3_dual_comparison.png)
Figure 4: 







##1.5 Conclusion


