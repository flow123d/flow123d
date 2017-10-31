#!/bin/bash
# Simple script which will execute input_convert.py file with python 3 interpret

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# populated by configure_file call
python3 ${SCRIPT_DIR}/yaml_converter/yaml_converter.py "$@"
