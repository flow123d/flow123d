#!/bin/bash

# Tar the build dir in order to use it in other jobs
# No compression is done as it seems to be slower. ~

build_dir=$1
cp config.cmake ${build_dir}/_config.cmake
tar -cvf build_dir.tar ${build_dir}/*
