#!/bin/sh
# 
# Copyright (C) 2007 Technical University of Liberec.  All rights reserved.
#
# Apply conversion on all YAML files store in subdirectories of '/tests/' directory,
#
# Conversion is performed by yaml_converter. See on submodule yaml_converter or on
# https://github.com/flow123d/yaml_converter

for dir in ../tests/*/     # list directories in the form "../tests/dirname/"
do
    dir=${dir%} 
    path=${dir##}/*.yaml
    echo "Conversion yaml files in directory ${path##}"
    python3 ./bin/yaml_converter/yaml_converter.py 4.0.0a01 ${path}
done
