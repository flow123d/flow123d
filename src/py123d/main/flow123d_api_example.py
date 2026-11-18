#!/bin/python3
# author: David Flanderka
#
# Example of Flow123d Python API using.


from dataclasses import dataclass
from py123d.main import Flow123dConfigBase


# Define Flow123d input source tree
input_data = {
    "flow123d_version": "4.0.0a01",

    "problem": {
        "TYPE": "Coupling_Sequential",
        "description": "Test8 - Steady flow with sources",
        "mesh": {
            "mesh_file": "./tests/00_mesh/square_1x1_shift.msh"
        },
        "flow_equation": { 
            "TYPE": "Flow_Darcy_LMH",
            "nonlinear_solver": {
                "linear_solver": {
                    "TYPE": "Petsc",
                    "r_tol": 1.0e-10,
                    "a_tol": 1.0e-10
                }
            },
            "input_fields": [
                { 
                    "region": "plane",
                    "anisotropy": 1,
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
            ],
            "output": {
                "fields": [
                    "pressure_p0",
                    "velocity_p0"
                ]
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


@dataclass
class MyFlowConfig(Flow123dConfigBase):
    """
    Flow123dConfigBase descendant, allows rewrite of environment arguments.
    Only 'output_dir' is changed, other arguments use parent values.
    """
    output_dir: str = "../test-result"


# Execute Flow123d simulation
config = MyFlowConfig()
config.run(input_data)
