#!/bin/python3
# author: David Flanderka.


from pathlib import Path
from dataclasses import dataclass, field
from py123d.main import Flow123dConfigBase


SCRIPT_DIR = Path(__file__).resolve().parent


# Define Flow123d input source tree
input_data = {
    "flow123d_version": "4.0.0a01",

    "problem": {
        "TYPE": "Coupling_Sequential",
        "description": "Test of FieldPython used in anisotropy and water_source_density (check results by switching to FieldFormula).",
        "mesh": {
            "mesh_file": str(SCRIPT_DIR / "../00_mesh/square_1x1_shift.msh")
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
                    "anisotropy": { 
                        "TYPE": "FieldPython",
                        "source_file": str(SCRIPT_DIR / "../02_generic_input/field_python_test.py"),
                        "class": "FieldPythonTest",
                        "used_fields": ["X"]

                    },
                    "water_source_density": { 
                        "TYPE": "FieldPython",
                        "source_file": str(SCRIPT_DIR / "../02_generic_input/field_python_test.py"),
                        "class": "FieldPythonTest",
                        "used_fields": ["X", "cross_section"]
                    },
                    "conductivity": { 
                        "TYPE": "FieldPython",
                        "source_file": str(SCRIPT_DIR / "../02_generic_input/field_python_test.py"),
                        "class": "FieldPythonTest",
                        "used_fields": ["X"]
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
class TestConfig(Flow123dConfigBase):
    """
    Flow123dConfigBase descendant, allows rewrite of environment arguments.
    """
    input_dir: str = field(
        default_factory=lambda: str(SCRIPT_DIR)
    )
    output_dir: str = field(
        default_factory=lambda: str(SCRIPT_DIR / "test_results/02_field_python.1")
    )

# Execute Flow123d simulation
config = TestConfig()
config.run(input_data)
