#!/bin/python3
# author: David Flanderka
#
# Example of Flow123d Python API using.


from pathlib import Path
from dataclasses import dataclass, field
from py123d.main import Flow123dConfigBase


SCRIPT_DIR = Path(__file__).resolve().parent


# Define Flow123d input source tree
input_data = {
    "flow123d_version": "4.0.0a01",

    "problem": {
        "TYPE": "Coupling_Sequential",
        "description": "Test8 - Steady flow with sources",
        "mesh": {
            "mesh_file": "./tests/00_mesh/square_1x1_shift.msh"                             # relative path to input_dir
            # "mesh_file": str(SCRIPT_DIR / "../../../tests/00_mesh/square_1x1_shift.msh")  # variant declaration of previous
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
class ExampleConfig(Flow123dConfigBase):
    """
    Flow123dConfigBase descendant, allows rewrite of environment arguments.
    Only 'input_dir', 'output_dir' are changed, other arguments use parent values.
    """
    input_dir: str = field(
        default_factory=lambda: str(SCRIPT_DIR / "../../..")
    )
    output_dir: str = field(
        default_factory=lambda: str(SCRIPT_DIR / "../../../test-result")
    )

# Execute Flow123d simulation
config = ExampleConfig()
config.run(input_data)
