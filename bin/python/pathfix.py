__author__ = 'Jan Hybs'

import sys, os

# print debug information about python
print ("Python:    " + str (sys.version).replace("\n", "") + ", " + str (sys.executable))

# add src/python into module path
path = os.path.join('..', '..', 'src', 'python')
sys.path.append (path)
