#!/usr/bin/env python3

import sys
import gmsh

gmsh.initialize(sys.argv, run=True)
gmsh.finalize()
