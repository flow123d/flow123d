#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
# 
# Simple python script testing embedded python functionality

def test():
    import os, sys, site, json, binascii
    # basic core modules
    modules = [os, site, json, binascii]
    # try to get __file__ property
    module_files = [getattr(x, '__file__', '') for x in modules]
    # return all existing __file__ properties joined by newline
    return '\n'.join([p for p in module_files if p])