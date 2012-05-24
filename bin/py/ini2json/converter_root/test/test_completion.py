# -*- coding: utf-8 -*-
'''
Created on 8.2.2012

@author: Tomáš Košek
'''
import unittest
import completion as comp

class Test(unittest.TestCase):

    def test_array_check(self):
        text = '{\n  "key": {\n    "1": {\n    },\n    "2": {\n    }\n  }\n}'
        res = '{\n  "key": [\n    {\n    },\n    {\n    }\n  ]\n}'
        self.assertEquals(res, comp.array_check(text))
        
    def test_value_type1(self):
        val  = '123.45'
        res = 123.45
        val_type = 'double'
        self.assertEquals(res, comp.value_type(val, val_type))
        
    def test_value_type2(self):
        val  = '1.5'
        res = 1
        val_type = 'int'
        self.assertEquals(res, comp.value_type(val, val_type))
    
    def test_value_type3(self):
        val  = '2'
        res = 'word2'
        val_type = 'enum(1 > word1, 2 > word2, 3 > word3)'
        self.assertEquals(res, comp.value_type(val, val_type))
        

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()