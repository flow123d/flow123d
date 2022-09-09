import flowpy

class PyFieldAnisotropy(flowpy.PythonFieldBase):
    pass

    def used_fields(self):
        field_list = ["X"]
        return field_list
    
    def __call__(self):
        diag = (0.5+2*z if (z < 0.25) else 1)
        return ( diag, 0, 0, 0, diag, 0, 0, 0, diag )


class PyFieldWaterSourceDensity(flowpy.PythonFieldBase):
    pass

    def used_fields(self):
        field_list = ["X"]
        return field_list
    
    def __call__(self):
        x = X[0]
        y = X[1]
        return ( 2*(1-x**2)+2*(1-y**2), )
