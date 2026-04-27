# encoding: utf-8
# author:   Jan Brezina

"""
This script check for older version of the IST output format and convert
it into the new format.

Syntax:
    python normalize_ist.py --version <version> <ist_file>

Performed actions:
- substitute "type_name" -> "name" (issue of old Record)
- "AbstractRecord" -> "Abstract"
- convert default_descentent value from the record name to its hash
- add wrapping header with the version info, the version could be passed by <version> parameter


"""

def sort_old(file):
    with open(file) as data_file:    
        data = json.load(data_file)       
        data.sort(key=sort_key)
        
        with open("sorted_"+file, "w") as output_file:
            json.dump(data, output_file, sort_keys=True, indent=4, separators=(',', ': '))
        
def sort_new(file):
    with open(file) as data_file:    
        data = json.load(data_file)
        node_list = data['ist_nodes']
        node_list.sort(key=sort_key)
        
        with open("sorted_"+file, "w") as output_file:
            json.dump(node_list, output_file, sort_keys=True, indent=4, separators=(',', ': '))