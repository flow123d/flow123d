import flowpy

class PythonAsm(flowpy.PythonFieldBase):
    pass

    # Method allows to define list of fields used in evaluation. This is the default defintion
    # that returns coords field and should be overwrite in descendant class. 
    def used_fields():
        field_list = ["velocity"]
        return field_list
    
    def vec_func(self):
        print("Calling call")
        return self.velocity * 2
