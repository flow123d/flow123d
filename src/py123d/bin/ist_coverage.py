#!/usr/bin/python
#
# Usage:
#   python ist_coverage.py <ist_file>.json
#
# This script is supposed to be run from the root dir.
# A single parameter is the IST description file.
#
# The script pass through all integration test yaml files and for every node in the IST tree 
# report files where the node is used. Nodes with different address are different so every type may appear there multiple times.
#
# Output is written into:
output_file = "ist_coverage_report.txt"

import yaml
import json
import sys
import os

class Instance:
    def __init__(self, tag, value):
        self.tag=tag
        self.value=value
        

def construct_record_instance_yaml(tag, loader, node):
    if isinstance(node, yaml.ScalarNode):
        value=loader.construct_scalar(node)
    if isinstance(node, yaml.SequenceNode):
        value=loader.construct_sequence(node)
    if isinstance(node, yaml.MappingNode):
        value=loader.construct_mapping(node)
    tag=str(tag)
    while tag[0] == "!":
        tag = tag[1:]
    return Instance(tag, value) 
    

class IST:
    @staticmethod
    def ist_node_by_id(node_id):
        return IST.type_dir[node_id]
    
    @staticmethod    
    def make_type(ist_node):        
        assert 'input_address' in ist_node, ist_node
        assert len(ist_node['input_address'])< 30, ist_node['input_address']
        input_type=ist_node['input_type']
        if input_type in ['String', 'Integer', 'Double', 'Selection', 'Bool', 'FileName']:
            input_type="Scalar"
        type_class=getattr(IST, input_type)        
        type_inst = type_class(ist_node)
        return type_inst
    
    
    class BaseType:
        def __init__(self, ist_type):
            self.ist=ist_type
            self.files={}
            assert 'input_address' in self.ist, self.ist
            assert self.ist['input_address'] != None, self.ist
            pass
        
        def append_file(self, filename):
            if filename in self.files:
                self.files[filename]+=1
            else:
                self.files[filename]=0
        
        def dfs_report_files_base(self, of):
            address = "/" + "/".join([ str(s).encode('utf-8') for s in self.ist['input_address'] ])
            files=list(self.files.keys())
            of.write("{}; {}; {}\n".format(address, len(files), files) )    
            
        def node_with_address(self, node, key):
            new_node=node.copy()
            new_node['input_address'] = self.ist['input_address'] + [ key ]
            return new_node
        
        def is_none(self, yaml):
            if (yaml == None):
                return False
            if (isinstance(yaml, str) and yaml.strip() == ""  ):
                return False
            self.yaml=yaml
            return True
            
        def consume_instance(self, yaml):
            '''
            If yaml is Instance object, set self.tag and return yaml.value.
            Otherwise return (None, yaml) pair.
            Also check for none values and possibly return None.
            '''
            if isinstance(yaml, Instance):
                self.tag = yaml.tag
                return self.is_none(yaml.value)
            else:
                self.tag=None
                return self.is_none(yaml)
    '''
    TODO:
    - call consume yaml at begin of add_file and transpose methods
    - in auto conversions pass down Instance(tag, yaml) instead of yaml itself.
    - try to meerge add_file and transpose a bit 
    '''
            
    class Scalar(BaseType):
        def __init__(self, ist_type):
            IST.BaseType.__init__(self, ist_type)
           
        def get_childs():
            return []
        
        def add_file(self, yaml, filename):
            if not self.consume_instance(yaml):
                return
            else:
                self.append_file(filename)
    
        def dfs_report_files(self, of):
            self.dfs_report_files_base(of)
            
        def transpose(self, yaml):
            if not self.consume_instance(yaml):
                return
            if type(self.yaml)==list:
                return self.yaml
            else:    
                return [ self.yaml ]
        
    class Array(BaseType):
        def __init__(self, ist_type):
            IST.BaseType.__init__(self, ist_type)
            node=IST.ist_node_by_id(self.ist['subtype'])
            new_node=self.node_with_address(node, 0)
            self.childs = [ IST.make_type(new_node) ]
            
        def add_file(self, yaml, filename):
            if not self.consume_instance(yaml):
                return
            self.append_file(filename)
            if type(self.yaml) != list:
                # autoconversion
                self.yaml=self.childs[0].transpose(Instance(self.tag,self.yaml))
            for item in self.yaml:
                self.childs[0].add_file(item, filename)
                    
            
        def dfs_report_files(self, of):
            self.dfs_report_files_base(of)
            self.childs[0].dfs_report_files(of)
        
        def transpose(self, yaml):
            if not self.consume_instance(yaml):
                return

            if type(self.yaml) == list:
                new_array=[]
                max_len=0
                for item in self.yaml:
                    t_item=self.childs[0].transpose(item)
                    new_array.append(t_item)
                    max_len=max(max_len, len(t_item))
                t_array=[]
                for i in range(max_len):
                    tt=[]
                    for t_item in new_array:
                        if len(t_item) > 1:
                            assert i<len(t_item), "Different lengths of arrays to transpose."
                            tt.append(t_item[i])
                        else:
                            tt.append(t_item[0])
                    t_array.append(tt)
                    
                return t_array
            else:    
                return [ [item] for item in self.childs[0].transpose(Instance(self.tag, self.yaml)) ]
                
            
    class Record(BaseType):
        def __init__(self, ist_type):
            IST.BaseType.__init__(self, ist_type)
            self.register_to_yaml()
            self.childs={}
            for key in self.ist['keys']:
                key_name=key['key']
                if key_name == "TYPE":
                    continue
                node=IST.ist_node_by_id(key['type'])
                new_node=self.node_with_address(node, key_name)
                self.childs[key_name]=IST.make_type(new_node)
                                 
        
        def register_to_yaml(self):
            tag="!"+self.ist['name']
            def constructor(loader, node):
                return construct_record_instance_yaml(tag, loader, node)            
            yaml.add_constructor(tag, constructor)

        def add_file(self, yaml, filename):
            if not self.consume_instance(yaml):
                return
                
            self.append_file(filename)
            if (type(self.yaml) == dict): 
                for (key, child) in self.childs.items():
                    if key in self.yaml:                        
                        child.add_file(self.yaml[key], filename)
                        
            else:
                # autoconversion
                key = self.ist['reducible_to_key']
                assert key and (key in self.childs), "Error: Invalid autoconversion key: {} in Record {}.".format(key, self.ist['name'])
                self.childs[key].add_file(Instance(self.tag, self.yaml), filename)

        def dfs_report_files(self, of):
            self.dfs_report_files_base(of)
            for (key, child) in self.childs.items():
                child.dfs_report_files(of)
                    
        def transpose(self, yaml):
            if not self.consume_instance(yaml):
                return

            if (type(self.yaml) == dict):                 
                new_dict={}
                max_len=0
                for (key, value) in self.yaml.items():
                    child=self.childs[key]
                    new_dict[key]=child.transpose(value)
                    max_len=max(max_len, len(new_dict[key]))
                
                t_dict=[]
                for i in range(max_len):
                    tt={}
                    for (key, t_item) in new_dict.items():
                        if len(t_item) > 1:
                            assert i<len(t_item), "Different lengths of arrays to transpose."
                            tt[key]=t_item[i]
                        else:
                            tt[key]=t_item[0]
                    t_dict.append(tt)
                return t_dict
            else:
                # autoconversion (need to be fixed as in ..., somehow merge these two passes)
                key = self.ist['reducible_to_key']
                
                assert key and (key in self.childs), "Error: Invalid autoconversion key: {} in Record {}.".format(key, self.ist['name'])
                return [ { key: t_item } for t_item in self.childs[key].transpose(Instance(self.tag, self.yaml)) ]

    class Tuple(Record):
        def __init__(self, ist_type):
            IST.Record.__init__(self, ist_type)
            self.child_list=[]
            for key in self.ist['keys']:
                key_name=key['key']
                if key_name == "TYPE":
                    continue
                node=IST.ist_node_by_id(key['type'])
                new_node=self.node_with_address(node, key_name)
                self.child_list.append( (key_name, IST.make_type(new_node) ) )
            
        def add_file(self, yaml, filename):
            if not self.consume_instance(yaml):
                return
                            
            if (type(self.yaml) == list): 
                if len(self.child_list)==len(self.yaml):
                    yaml_dict={}
                    for (key_type, value) in zip(self.child_list, self.yaml):
                        yaml_dict[ key_type[0] ] = value
                    self.yaml=yaml_dict    
                else:
                    assert false, "Defaults for incomplete tuples not implemented yet."
            IST.Record.add_file(self, Instance(self.tag, self.yaml), filename)

        def transpose(self, yaml):
            if not self.consume_instance(yaml):
                return
                            
            if (type(self.yaml) == list): 
                if len(self.child_list)==len(self.yaml):
                    yaml_dict={}
                    for (key_type, value) in zip(self.child_list, self.yaml):
                        yaml_dict[ key_type[0] ] = value
                    self.yaml=yaml_dict    
                else:
                    assert false, "Defaults for incomplete tuples not implemented yet."
            return IST.Record.transpose(self, Instance(self.tag, self.yaml))

    class Abstract(BaseType):
        def __init__(self, ist_type):
            IST.BaseType.__init__(self, ist_type)
            self.childs={}
            for rec in self.ist['implementations']:
                rec_node=IST.ist_node_by_id(rec)
                rec_name=rec_node['name']                
                new_node=self.node_with_address(rec_node, rec_name)
                self.childs[rec_name]=IST.make_type(new_node)                
                
                
        def add_file(self, yaml, filename):
            if not self.consume_instance(yaml):
                return

            self.append_file(filename)
                       
            if self.tag:
                assert self.tag in self.childs.keys(), "Error: Wrong tag {} for abstract {} at address {}.".format(self.tag, self.ist['name'], self.ist['input_address'])
                self.childs[self.tag].add_file(self.yaml, filename)
            else:
                # autoconversion
                try:
                    self.tag=IST.ist_node_by_id( self.ist['default_descendant'] )['name']
                except:
                    assert False, "Missing tag for Abstract: {} at address {}.".format(self.ist['name'], self.ist['input_address'])
                    
                self.childs[self.tag].add_file(self.yaml, filename)


        def dfs_report_files(self, of):
            self.dfs_report_files_base(of)
            for child in self.childs.values():
                child.dfs_report_files(of)
        
        def transpose(self, yaml):
            if not self.consume_instance(yaml):
                return
            if not self.tag:
                try:
                    self.tag=IST.ist_node_by_id( self.ist['default_descendant'] )['name']
                except:
                    assert False, "Missing tag for Abstract: {} at address {}.".format(self.ist['name'], self.ist['input_address'])
            
            assert self.tag in self.childs.keys(), "Error: Wrong tag {} for abstract {} at address {}.".format(self.tag, self.ist['name'], self.ist['input_address'])
            self.yaml = self.childs[self.tag].transpose(self.yaml)
            return [ Instance(self.tag, item) for item in self.yaml ]

    
def make_ist(ist_flat):
    types=ist_flat['ist_nodes']
    IST.type_dir={}
    
    for ist_type in types:
        # skip last empty dict
        if not ist_type:
            continue
        IST.type_dir[ist_type['id']]=ist_type
        if ist_type['name'] == 'Root':
            root_type=ist_type
    
    root_type['input_address']=[]
    return IST.Record(root_type)
        
    
####################### MAIN


ist_file = sys.argv[1]
with open(ist_file,"r") as f:
    ist_flat=json.load(f)

root_type = make_ist(ist_flat)

walk_root_dir = "./tests"
count_files=0
for subdir, dirs, files in os.walk(walk_root_dir):
    for f in files:
        if ( f.endswith('yaml')
            and f != "config.yaml"
            and not "ref_out" in str(subdir)
            and not "test_results" in str(subdir) ):                
                with open(os.path.join(subdir, f)) as ff:
                    raw_yaml=yaml.load(ff)
                count_files+=1
                filename=os.path.join(os.path.relpath(subdir, walk_root_dir), f)                
                print( "Processing: " + filename )
                root_type.add_file(raw_yaml, filename)

print("Total tests: ", count_files)          
with open(output_file, "wt") as of:
    root_type.dfs_report_files(of)
print("Coverage report written into: " + output_file)
