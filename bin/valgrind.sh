#!/bin/bash

# call flow123d under valgrind
valgrind --suppressions=${0%/*}/python.supp ${0%/*}/flow123d "$@"