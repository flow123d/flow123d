__author__ = 'Jan Hybs'
import sys, os



def require_version_2 ():
    if (sys.version_info > (3, 0)):
        print ('Error: Python 2 is required')
        sys.exit (1)

def require_version_3 ():
    if (sys.version_info < (3, 0)):
        print ('Error: Python 3 is required')
        sys.exit (1)
