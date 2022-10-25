import flowpy
import numpy as np

class FieldPythonTest(flowpy.PythonFieldBase):
    pass

    def water_source_density(self):
        XY = self.X[0:2]
        AB = 2*(1-XY**2)
        return AB[0] + AB[1] + self.cross_section   
            # or np.sum(AB, axis=0) + self.cross_section

    def anisotropy(self):
        z = self.X[2]
        diag = np.where( z < 0.25, 0.5+2*z, 1)
        return self.repl(np.eye(3)) * diag
