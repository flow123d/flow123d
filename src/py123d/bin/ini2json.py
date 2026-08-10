#!/usr/bin/python

'''
This script is like a link. It calls converter from *.ini input
into modified JSON. For syntax call:

python ini2json.py --help

Conversion of the INI file is described by the template file:
"ini2json/converter_root/converter/ini_to_json_template.txt"

The program can also merge a *.mtr file into created JSON file.
However, this is hardwired in 
"ini2json/converter_root/converter/materials.py"

'''

import sys, os

conv_path= sys.path[0] + os.sep +"ini2json" + os.sep + "converter_root" + os.sep + "converter"
print conv_path

sys.path.append(conv_path)
import converter 

converter.main()