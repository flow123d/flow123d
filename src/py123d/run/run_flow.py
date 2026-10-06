#!/bin/python3
# author: David Flanderka

import flow123d

input_data = {
    "flow123d_version": "4.0.0a01",

    "problem": {
        "TYPE": "Coupling_Sequential",
        "description": "Test8 - Steady flow with sources",
        "mesh": {
            "mesh_file": "../00_mesh/square_1x1_shift.msh"
        }.
        "flow_equation": { 
            "TYPE": "Flow_Darcy_LMH",
            "nonlinear_solver": {
                "linear_solver": {
                    "TYPE": "Petsc",
                    "r_tol": 1.0e-10,
                    "a_tol": 1.0e-10
                }
            },
            "input_fields": {
                { 
                    "region": "plane",
                    "anisotropy": 1
                    "water_source_density": { 
                        "TYPE": "FieldFormula",
                        "value": "2*(1-X[0]**2)+2*(1-X[1]**2)"
                    }
                 },
                 {
                     "region": ".plane_boundary",
                     "bc_type": "dirichlet",
                     "bc_pressure": 0
                 }
            },
            "output": {
                "fields": {
                    "pressure_p0",
                    "velocity_p0"
                }
            },
            "output_stream": {
                "file": "./flow.pvd",
                "format": {
                    "TYPE": "vtk",
                    "variant": "ascii"
                }
            }
        }
    }
}

flow123d.run(input_data)
