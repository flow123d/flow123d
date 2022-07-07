#!/bin/bash
#
# This script adds paths to pybind repository and to flowpy library to PYTHONPATH variable and runs Python.

# Path to pybind
export PYTHONPATH="${PYTHONPATH}:third_party/pybind11-2.9.2/include"
# Path to actual building directory
export PYTHONPATH="${PYTHONPATH}:build-DF_asm_field_python/src"
# Set actual Python instalation
python3.8
