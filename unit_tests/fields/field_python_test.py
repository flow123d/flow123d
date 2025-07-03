from py123d import field as flowpy
import numpy as np

class PythonAsm(flowpy.PythonFieldBase):
    pass

    def vec_func(self):
        print("Calling call")
        return self.velocity * 2


class FieldPythonTest1(flowpy.PythonFieldBase):
    pass

    def scalar_field(self):
        """ Evaluates expression: x """
        return self.X[0]

    def vector_field(self):
        """ Evaluates expression: [x, 2*x, 0.5] """
        xx = self.X[0]
        part_xy = self.repl( np.array([1.0, 2.0, 0.0]) ) * xx
        part_z  = self.repl( np.array([0.0, 0.0, 0.5]) ) * np.full((self._region_chunk_end-self._region_chunk_begin), 1.0)
        return part_xy + part_z

    def tensor_field(self):
        """ Evaluates expression: [[x, 0.2, 0.3], [0.2, 0.4, 0.5], [0.3, 0.5, 0.6]] """
        xx = self.X[0]
        part_00 = self.repl( np.array([[1.0, 0, 0], [0, 0, 0], [0, 0, 0]]) ) * xx
        other   = self.repl( np.array([[0.0, 0.2, 0.3], [0.2, 0.4, 0.5], [0.3, 0.5, 0.6]]) ) * np.full((self._region_chunk_end-self._region_chunk_begin), 1.0)
        return part_00 + other


class FieldPythonTest2(flowpy.PythonFieldBase):
    pass

    def scalar_field(self):
        """ Evaluates expression: y """
        return self.X[1]

    def vector_field(self):
        """ Evaluates expression: [y, 2*y, 0.5] """
        yy = self.X[1]
        part_xy = self.repl( np.array([1.0, 2.0, 0.0]) ) * yy
        part_z  = self.repl( np.array([0.0, 0.0, 0.5]) ) * np.full((self._region_chunk_end-self._region_chunk_begin), 1.0)
        return part_xy + part_z

    def tensor_field(self):
        """ Evaluates expression: [[y, 2.2, 2.3], [2.2, 2.4, 2.5], [2.3, 2.5, 2.6]] """
        yy = self.X[1]
        part_00 = self.repl( np.array([[1.0, 0, 0], [0, 0, 0], [0, 0, 0]]) ) * yy
        other   = self.repl( np.array([[0.0, 2.2, 2.3], [2.2, 2.4, 2.5], [2.3, 2.5, 2.6]]) ) * np.full((self._region_chunk_end-self._region_chunk_begin), 1.0)
        return part_00 + other


class FieldPythonTest3(flowpy.PythonFieldBase):
    pass

    def python_field(self):
        """ Evaluates expression: x """
        return self.X[0]
        