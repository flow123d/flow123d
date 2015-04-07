# -*- coding: utf-8 -*-
'''
Created on 8.2.2012

@author: Tomáš Košek
'''
import unittest
import conversion as cnv

class Test(unittest.TestCase):


    def test_change_words(self):
        word = "on"
        res = "true"
        self.assertEquals(res, cnv.change_words(word))
        
        
    def test_which_type(self):
        value = "1e-11"
        res = 1
        self.assertEquals(res, cnv.which_type(value))
        
        
    def test_read_header(self):
        value = "[header]"
        res = '\t"header" : {\n'
        self.assertEquals(res, cnv.read_header(value))


    def test_read_data(self):
        value = "value = string" #"value = 10"
        res = '"value" : "string"' #'"value" : 10'
        self.assertEquals(res, cnv.read_data(value))


    def test_read_array(self):
        value = "A B C 1 2 3"
        res = '["A", "B", "C", 1, 2, 3]'
        self.assertEquals(res, cnv.read_array(value))

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()