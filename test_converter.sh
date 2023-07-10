#!/bin/sh

for dir in ./tests/*/     # list directories in the form "./tests/dirname/"
do
    dir=${dir%} 
    path=${dir##}/*.yaml
    echo "Conversion yaml files in directory ${path##}"
    python3 ./bin/yaml_converter/yaml_converter.py 4.0.0 ${path}
done
