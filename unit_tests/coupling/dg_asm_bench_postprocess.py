"""
Script selects assembly data from profiler output and converts its to csv table.
Usage:

    python3 profiler_to_csv.py <profiler_file>
"""

import sys
import json
import csv
import pandas as pd


"""
Hold processed profiler data and helper variables. 
Object is used during processing of profiler JSON input.
"""
class ProfilerHandler:
    def __init__(self, profiler_file, program_params):
        self.profiler_file = profiler_file
        self.program_params = program_params
        
        # Define name of columns 
        self.column_names = ['commit', 'run_id', 'domain_shape', 'mesh_size', 'spacedim', 'uniformity', 'n_mesh_elements', 'n_repeats', 'field_variant',   
            'assembly_variant', 'assembly_class', 'tag', 'mumber_of_calls', 'integral_type', 'time', 'time_fraction', 'time_fraction_of_reminder']
        
        # Dictionary holds pairs tag-method, these methods are called for selected tags during processing of profiler tree
        self.tag_method = {
            'full_mesh':                     'full_mesh', 
            'BaseMeshReader - mesh factory': 'empty', 
            'ZERO-TIME STEP':                'empty', 
            'MassAssembly':                  'assembly', 
            'StiffnessAssembly':             'assembly', 
            'SourcesAssembly':               'assembly'} 
        
        # Create other data members using during profiler processing
        self.json_tree_nodes = []
        self.n_mesh_elements = 0
        self.mesh_repeats = 0
        self.common_row_data = []
        self.sum_subtags_time = 0.0


# Process one assembly (Mass, Stiffness or Sources) and all subtags 
def assembly(json_node, ph, df):
    integrals = {
        "assemble_volume_integrals": "cell",
        "assemble_fluxes_boundary":  "boundary_side",
        "assemble_fluxes_elem_elem": "edge",
        "assemble_fluxes_elem_side": "dimjoin"
    }

    ph.json_tree_nodes.append(json_node)
    asm_name = json_node['tag']
    #print(" -", json_node['tag'], ph.json_tree_nodes[-2]['tag'], ph.json_tree_nodes[-1]['tag'])
    ph.sum_subtags_time = 0.0
    for i in json_node['children']:
        tag = i['tag']
        time = i['cumul-time-sum']
        integral = None
        if tag in integrals:
            integral = integrals[tag]
        tag_row_data = ph.common_row_data.copy()
        tag_row_data.extend( [asm_name, tag, i['call-count-sum'], integral, time, (time/json_node['cumul-time-sum']), 0] )
        df = pd.concat([df, pd.DataFrame([tag_row_data], columns=ph.column_names)], ignore_index=True)
        ph.sum_subtags_time += time

    asm_time = json_node['cumul-time-sum']
    asm_row_data = ph.common_row_data.copy()
    asm_row_data.extend( [asm_name, asm_name, json_node['call-count-sum'], None, asm_time, 0, ((asm_time-ph.sum_subtags_time)/asm_time)] )
    df2 = pd.concat([df, pd.DataFrame([asm_row_data], columns=ph.column_names)], ignore_index=True)
    ph.json_tree_nodes.pop()
    return df2


# We can skip some tags which are not to be processed (e.g. 'ZERO-TIME STEP')
def empty(json_node, ph, df):
    return df


# Process 'full_mesh' tag
def full_mesh(json_node, ph, df):
    ph.mesh_repeats = json_node['call-count-sum']
    ph.n_mesh_elements = json_node['children'][1]['call-count-sum']
    mesh_full_name = ph.json_tree_nodes[-3]['tag'] + '_' + ph.json_tree_nodes[-2]['tag'] + '_' + ph.json_tree_nodes[-1]['tag']
    mesh_params = mesh_full_name.split('_')
    ph.common_row_data = [ ph.program_params[0], ph.program_params[1], mesh_params[0], mesh_params[3], mesh_params[1], 
                    mesh_params[2], (ph.n_mesh_elements/ph.mesh_repeats), ph.mesh_repeats, mesh_params[5], mesh_params[4] ]
    mesh_row_data = ph.common_row_data.copy()
    mesh_row_data.extend( ['', 'full_mesh', json_node['call-count-sum']*ph.mesh_repeats, '',  json_node['cumul-time-sum'], '', ''] )
    df2 = pd.concat([df, pd.DataFrame([mesh_row_data], columns=ph.column_names)], ignore_index=True)
    return process_node(json_node, ph, df2)


def process_node(json_node, ph, df):
    ph.json_tree_nodes.append(json_node)
    if 'children' in json_node:
        for i in json_node['children']:
            if i['tag'] in ph.tag_method.keys():
                possibles = globals().copy()
                possibles.update(locals())
                method = possibles.get( ph.tag_method[ i['tag'] ] )
                df = method(i, ph, df) 
            else:
                df = process_node(i, ph, df)
    ph.json_tree_nodes.pop()
    return df


# Load data from profiler JSON file
def load_profiler_data(ph, df):
    # Load JSON file
    with open(ph.profiler_file) as f_in:
        profiler_data = json.load(f_in)
    
    # Only one descendant 'Whole Program' of root tree exists 
    whole_program = profiler_data['children'][0]
    
    df = process_node(whole_program, ph, df)
    return df


# Perform outout to CSV file
def csv_output(csv_file, program_params, df):
    out_mode = 'a'
    out_header = False
    if program_params[1] == '1': 
        out_mode = 'w'
        out_header = True
    df.to_csv(path_or_buf=csv_file, header=out_header, index=False, mode=out_mode)


def main():
    pass
    
    # store names of input json and output csv file from arg
    profiler_file = sys.argv[1]
    # commit, run_id 
    program_params = sys.argv[2:4]
    csv_file = profiler_file + '.csv'

    # load and process JSON
    ph = ProfilerHandler(profiler_file, program_params)
    df = pd.DataFrame( columns=ph.column_names )
    df = load_profiler_data(ph, df)
    
    # output to CSV
    csv_output(csv_file, program_params, df)

# Processes one node of profiler tree: print tag and calls children nodes recursivelly
def print_profiler_node(json_in, indent):
    print(indent, json_in['tag'])
    if 'children' in json_in:
        indent2 = '    ' + indent
        for i in json_in['children']:
            print_profiler_node(i, indent2)

"""
Helper function. Prints tree of profiler nodes.

Usage:
  python3 ./dg_asm_bench_postprocess.py <profiler_output.json> <commit> <run_id>
"""
def print_profiler_tree():
    profiler_file = sys.argv[1]
    # Load JSON file
    with open(profiler_file) as f_in:
        profiler_data = json.load(f_in)
    print_profiler_node(profiler_data['children'][0], '-') # Whole program tag


if __name__ == "__main__":
    #print_profiler_tree()
    main()

