#!/bin/python3
# author: David Flanderka

from dataclasses import dataclass, field, fields
from typing import Any

import flow123d_python_api


@dataclass
class Flow123dConfigBase:
    """
    Base configuration for running Flow123d from Python.

    Derive your own configuration class and override
    only the parameters that you need.
    """

    # Environment arguments
    log_filename: str = field(
        default="",
        metadata={"help": "Log file name; '' means default, '\n' means no logging."}
    )
    petsc_redirect_file: str = field(
        default="",
        metadata={"help": "Optional file name for output of PETSc parameters."}
    )
    signal_handler_off: bool = field(
        default=False,
        metadata={"help": "Turn off signal handling useful to debug with valgrind."}
    )
    use_profiler: bool = field(
        default=True,
        metadata={"help": "If true, we do output of profiling information."}
    )
    memory_monitoring: bool = field(
        default=False,
        metadata={"help": "If true, memory monitoring is switched on."}
    )
    profiler_path: str = field(
        default="",
        metadata={"help": "Location of the profiler report file."}
    )
    input_dir: str = field(
        default=".",
        metadata={"help": "Input (working) directory."}
    )
    output_dir: str = field(
        default="./test-result",
        metadata={"help": "Output directory."}
    )
    input_filename: str = field(
        default="",
        metadata={"help": "If there is 'solve' option."}
    )

    # Interface to Flow123d

    def get_env_args(self) -> dict[str, Any]:
        """Convert configuration parameters to a Python dictionary."""

        return {
            item.name: getattr(self, item.name)
            for item in fields(self)
        }

    def run(self, input_data: dict):
        """Run Flow123d using the supplied input dictionary."""

        env_args = self.get_env_args()

        flow123d_python_api.run(
            input=input_data,
            env_args=env_args
        )

    @classmethod
    def print_help(cls):
        """Print available Flow123d environment parameters."""

        print("Flow123d environment parameters:\n")

        for item in fields(cls):
            help_text = item.metadata.get("help", "")
            print(f"  {item.name}")
            print(f"      default: {item.default}")
            print(f"      type:    {item.type}")
            print(f"      help:    {help_text}")
            print()
