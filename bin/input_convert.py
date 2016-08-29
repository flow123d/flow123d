#!/usr/bin/python3

'''
  Usage:    
    input_convert.py <files to convert>
    
  The script converts the input CON file for the version 1.8.2 or the YAML file for the version 2.0.0_rc
  the the input for the version 2.0.0. The original YAML files are stored as *.N.orig.yaml, where N is an
  increasing number avoiding to overwrite previous saved files.   
'''  


import sys
import os

# add common dir to python path
common_dir = os.path.join(os.path.split(
    os.path.dirname(os.path.realpath(__file__)))[0], "src/python/GeoMop/common")
sys.path.insert(1, common_dir)
# add importer dir to python path
importer_dir = os.path.join(os.path.split(
    os.path.dirname(os.path.realpath(__file__)))[0], "src/python/GeoMop/ModelEditor")
sys.path.insert(1, importer_dir)

import importer


def convert_file(in_file):
    if not os.path.exists(in_file):
      raise Exception("Can not find file: {} ".format(in_file))
    (root, ext) = os.path.splitext(in_file)
    
    con_file=""
    yaml_file=""
    if ext == ".con":
        input_version="1.8.2"
        con_file=in_file        
    elif ext == ".yaml":
        input_version="2.0.0_rc"
        # save original     
        N=0
        while (True):
            save_name = root + ".{}.orig.yaml".format(N)
            if not os.path.exists(save_name) :
                break
            N+=1
        
        # simple file copy
        with open(in_file, "r") as f_in:
            with open(save_name, "w") as f_out:
                for line in f_in:
                    f_out.write(line)
                    
        yaml_file = save_name            
    else:
        raise Exception("Unknown input file type: " + ext)
            
    dest_file= root + ".yaml"    
    if (os.path.exists(dest_file)):    
        os.remove(dest_file)
    
    output_version="2.0.0"    
    transform_name="from_"+input_version+"_to_"+output_version;                
    
    importer.do_transform(con_file, yaml_file, dest_file, [ transform_name ])
    
    
for in_file in sys.argv[1:] :
    convert_file(in_file)
    
  