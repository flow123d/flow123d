""" Usage example of FieldPython. 
    Define for correct functionality:
     - class of arbitrary name, this class must be descendant of flowpy.PythonFieldBase
     - method of same name as name of evaluated field
       - shape of return value must correspond to shape of result
     - one class can evaluate one or more fields
    Define following keys in declaration of FieldPython in YAML input file:
     - source_file: name of python file with extension
     - class: name of class 
     - used_fields: list of fields necessary for evaluation (see notes below) 
    
    author: David Flanderka
"""

import flowpy
import numpy as np

class FieldPythonTest(flowpy.PythonFieldBase):
    pass

    def water_source_density(self):
        """ Evaluates expression: 2*(1-x^2)+2*(1-y^2)+cross_section
            List of used_fields is: ["X", "cross_section"]
            (the value of X corresponds to coordinates) 
        """
        XY = self.X[0:2]
        AB = 2*(1-XY**2)
        return AB[0] + AB[1] + self.cross_section   
            # or np.sum(AB, axis=0) + self.cross_section

    def anisotropy(self):
        """ Evaluates expression: if(z < 0.25, 0.5+2*z, 1)
            (predefined method 'repl' ensures that resulting value is stored in all evaluated points) 
            List of used_fields is: ["X"]
        """
        z = self.X[2]
        diag = np.where( z < 0.25, 0.5+2*z, 1)
        return self.repl(np.eye(3)) * diag

    def conductivity(self):
        """ Evaluates expression: if(z < 0.25, 0.5+2*z, 1)
            List of used_fields is: ["X"]
        """
        A = np.array([[0.2, 0, 0], [0, 0.3, 0], [0, 0, 0.5]])
        return np.linalg.norm(A @ self.X, axis = 0)
