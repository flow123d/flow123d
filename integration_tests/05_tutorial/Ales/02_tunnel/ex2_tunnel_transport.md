
# FLOW123D EXAMPLES
## 1. Example 2 – tunnel
Again the example is inspirited by the water treatment plant tunnel Bedřichov in the granite rock massif and the particular seepage site 23 m under the surface. This seepage site has fast reaction of precipitation and various chemical values are measured.

### 1.1    Geometry and boundary conditions
The geometry was considered as rectangle 500 × 300 m with a octagon hole 23 meters under surface which represents the tunnel. So the tunnel was consider as infinite tube. 

    mesh:
      mesh_file: ./input/mesh.gmsh

## 1.2    Hydraulic model
The hydraulic model was fitted on the shape of flux field, where it was assumed that the tunnel drains only part of the model surface. And we fit the model on the estimated discharge of the seepage site.

### 1.2.2  Boundary condition

Thickness of model was set for 1 meter ("cross_section"). The hydraulic conductivity of rock medium was prescribed on 1e-7.
On the surface (".surface") the annual precipitation was prescribed (2.698E-08 m/s = 852 mm/yr), the pressure 270 m was prescribed on the base (".base") because of assumption of local groundwater flow. On the tunnel (".tunnel") the atmospheric pressure was set. Others boundaries was prescribed no flow automatically (Figure 1).  Code of control file "ex02_tunnel_transport.yaml" is:

    input_fields:
      - region: rock
        conductivity: 1.E-07
        cross_section: 1
      - region: .tunnel
        bc_type: dirichlet
        bc_pressure: 0
      - region: .base
        bc_type: dirichlet
        bc_pressure: 270
      - region: .surface
        bc_type: total_flux
        bc_flux: 6.34E-09
   

![alt text](http://bacula.nti.tul.cz/~ales.balvin/ex2_tunnel/ex2_bcd)

Figure 1: Geometry and boundary condition of model

### 1.2.2  Results
The results are depicted in Figures 1, where the flux field is shown with pressure (range from zero till max. value). The unsaturated layer is shown with peizometric field (file "flow.msh" postprocesed in GMSH).
The outflow through boundary ".tunnel" is 7.88048e-06 m/s (in file "water_balance.txt"). It is 7.8E-03 l/s. The mean discharge of seepage site is approximately similar 0.036 l/s/m when it is considered that the measured discharge is only 10% per cent of overall. 

![alt text](http://bacula.nti.tul.cz/~ales.balvin/ex2_tunnel/ex2_results_of_flux.png)

Figure 2: Pressure, boundary of water level and piezometric head in unsaturated zone and flux field. 

## 1.3   Transport model - Pulse
The transport model is inspirited by the measurement of stable isotopes O-18. From measurement the time of movement of water molecule from surface into the tunnel is estimated to 38 months. The porosity and other parameters for advection-dispersion equation was fitted on this time.

### 1.3.1  Model Settings
The settings of numerical diffusion ("transport: !Solute_Advection_FV") is similar like for Example 1. Except that the relative concertation 100% is prescribed only for first time step. The porosity was set in the some section. Then the code of control file "ex02_tunnel_transport.yaml" is:
   
    transport: !Solute_Advection_FV
      input_fields:
        - region: rock
          porosity: 0.067
        - region: .surface
          bc_conc: 100
          time: 0
        - region: .surface
          bc_conc: 0
          time: 1e7

The time step was set on 1e7 second and time end of simulation was set on 1e10.       
      
## 1.3.2  Results of Pulse transport
The output tunnel breakthrough curve is depicted on the Figure 3 and the movement of mass is depicted on the Figure 4. The transit time could be around 40 months, which is approximately the centre of gravity of the curve.

![alt text](http://bacula.nti.tul.cz/~ales.balvin/ex2_tunnel/ex2_results_of_tranpost_curve.png)

Figure 3: Breakthrough curve of mass at the tunnel. 

![alt text](http://bacula.nti.tul.cz/~ales.balvin/ex2_tunnel/ex2_results_of_transport.png)
 
Figure 4: The mass movement in the model rock massive.


## 1.4  Variant - transport of real isotopes
On the locality is sampled the stable isotope O-18 in monthly steps in precipitation at nearby experimental catchment Uhlirska and at the seepage site 23m depth. Previous breakthrough curve variant is used for computation of transient stable isotopes modelling.

### 1.4.1  Model Settings of transport of real isotopes
In the file "ex02_tunnel_transport_transient_transport.yaml" units in flow boundary condition are changed on m/days. On the ".tunnel" is prescribed outflow "bc_flux: -9.16E-02" (m/day).
Transport was computed in two numerical method. Fist is in finite volume method only with numerical diffusion ("transport: !Solute_Advection_FV") and then discontinuous Galerkin method ("transport: !Solute_AdvectionDiffusion_DG") with possibility prescribed diffusion and longitudinally and transversally dispersion.
Boths methods are spread by initial condition in "input_fileds". For transport a time series is used ("inC.txt"), where column is with dates of measurements (from 1/2006 till 6/2013) and the second column is concentration of measurement of stable isotopes O-18. The concentration is prescribed on the ".surface" in every 15th day in a month. The generator of this time series for yaml is prepared in excel ("inC.xlsx").


    transport: !Solute_AdvectionDiffusion_DG
      input_fields:
        - region: rock
          porosity: 0.067
          init_conc: -10.5
          diff_m: 6e-10
          disp_l: 5
          disp_t: 1
        - region: .surface
          bc_conc: -12.85443
          time: 11.8333333330665
        - region: .surface
          bc_conc: -14.00255
          time: 42.2499999997308
        - region: .surface
          bc_conc: -12.80081
          time: 72.666666666395
        - region: .surface
          bc_conc: -12.34748
          time: 103.083333333059
		...
		...


    transport: !Solute_Advection_FV
      input_fields:
        - region: rock
          porosity: 0.067
          init_conc: -10.5
        - region: .surface
          bc_conc: -12.85443
          time: 11.8333333330665
        - region: .surface
          bc_conc: -14.00255
          time: 42.2499999997308
        - region: .surface
          bc_conc: -12.80081
          time: 72.666666666395
        - region: .surface
          bc_conc: -12.34748
          time: 103.083333333059
		...
		...

### 1.4.2  Results of transport of real isotopes

Results of both numerical method gave identical results. Fit of measured data is quit reasonable (Figure 5). The movement of transport is depicted in Figure 6.

![alt text](http://bacula.nti.tul.cz/~ales.balvin/ex2_tunnel/ex2results-graph_real_transport.png)
 
Figure 5: Result of real isotopes transport of stable isotopes. Concentration os isotopes on seepage site 23m under the surface.

![alt text](http://bacula.nti.tul.cz/~ales.balvin/ex2_tunnel/ex2_results_transport.png)
 
Figure 6: Result of real isotopes transport in two-dimensional model.


## 1.5   Conclusion
This example inspirited real locality and measurement shown the two-dimensional model with a tunnel in particular depth. The results was fitted and prepared for this example. The transport with real isotope was shown with two different numerical methods.
