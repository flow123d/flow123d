"""
Script selects assembly data from profiler output and converts its to csv table.
Usage:

    python3 profiler_to_csv.py <profiler_file>
"""

import sys
import os
import json
import csv
import pandas as pd


"""
Hold processed profiler data and helper variables. 
Object is used during processing of profiler JSON input.
"""
class ProfilerHandler:
    def __init__(self, profiler_file):
        self.profiler_file = profiler_file
        
        # Define name of columns 
        self.column_names = ['branch', 'commit', 'run_id', 'domain_shape', 'mesh_size', 'spacedim', 'uniformity', 'n_mesh_elements', 'n_repeats', 'field_variant',   
            'assembly_variant', 'assembly_class', 'tag', 'number_of_calls', 'integral_type', 'time', 'time_fraction', 'time_fraction_of_reminder']
        
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
        self.program_branch = ''
        self.program_revision = ''
        self.run_id = 0

    # Set base program parameters
    def set_program_params(self, program_branch, program_revision, run_id):
        self.program_branch = program_branch
        self.program_revision = program_revision
        self.run_id = run_id


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
        integral = ''
        if tag in integrals:
            integral = integrals[tag]
        tag_row_data = ph.common_row_data.copy()
        tag_row_data.extend( [asm_name, tag, i['call-count-sum'], integral, time, (time/json_node['cumul-time-sum']), ''] )
        df = pd.concat([df, pd.DataFrame([tag_row_data], columns=ph.column_names)], ignore_index=True)
        ph.sum_subtags_time += time

    asm_time = json_node['cumul-time-sum']
    asm_row_data = ph.common_row_data.copy()
    asm_row_data.extend( [asm_name, asm_name, json_node['call-count-sum'], '', asm_time, '', ((asm_time-ph.sum_subtags_time)/asm_time)] )
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
    ph.common_row_data = [ ph.program_branch, ph.program_revision, ph.run_id, mesh_params[0], mesh_params[3], mesh_params[1], 
                    mesh_params[2], (ph.n_mesh_elements/ph.mesh_repeats), ph.mesh_repeats, mesh_params[5], mesh_params[4] ]
    mesh_row_data = ph.common_row_data.copy()
    mesh_row_data.extend( ['', 'full_mesh', json_node['call-count-sum'], '',  json_node['cumul-time-sum'], '', ''] )
    df2 = pd.concat([df, pd.DataFrame([mesh_row_data], columns=ph.column_names)], ignore_index=True)
    return process_node(json_node, ph, df2)


# Process one tag and call recursivelly processing of children
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


# Unify values of DataFrame selected columns.
def unify_df_values(df):
    # unify 'cube' and 'square' to 'box'
    shape_map = {'square': 'B', 'cube': 'B', 'lshape': 'L'}
    df['domain_shape'] = [shape_map[s] for s in df['domain_shape']]
    
    # use short formats of mesh size: 'S', 'M' and 'L'
    size_map = {'small': 'S', 'medium': 'M', 'big': 'L'}
    df['mesh_size'] = [size_map[s] for s in df['mesh_size']]
    
    # use name of assembly class without 'Assembly' postfix
    asm_map = {'MassAssembly': 'Mass', 'StiffnessAssembly': 'Stiffness', 'SourcesAssembly': 'Sources', '': ''}
    df['assembly_class'] = [asm_map[s] for s in df['assembly_class']]
    
    return df


"""
Load data from set of profiler JSON files.

Name of files must be in format '<profiler_file>_<n>.json' where:
    <profiler_file> Same for all files in set, it is given by arg 'profiler_file'
    <n> Order of file in set 1,2,3,... Number of files is given by arg 'n_runs'
"""
def load_profiler_data(profiler_file, n_runs):

    # Create ProfilerHandler and DataFrame
    ph = ProfilerHandler(profiler_file)
    df = pd.DataFrame( columns=ph.column_names )
    
    for run_id in range(1, n_runs+1):
        file_name = ph.profiler_file + '_' + str(run_id) + '.json'

        # Load JSON file
        with open(file_name) as f_in:
            profiler_data = json.load(f_in)
    
        # Only one descendant 'Whole Program' of root tree exists 
        whole_program = profiler_data['children'][0]
        program_branch = profiler_data['program-branch']
        program_revision = profiler_data['program-revision']
        ph.set_program_params(program_branch, program_revision, run_id)
    
        df = process_node(whole_program, ph, df)
    
    df = unify_df_values(df)
    return df


# Perform outout to CSV file
def csv_output(csv_file, run_id, df):
    out_mode = 'a'
    out_header = False
    if run_id == '1': 
        out_mode = 'w'
        out_header = True
    df.to_csv(path_or_buf=csv_file, header=out_header, index=False, mode=out_mode)


def main():
    pass
    
    # store names of input json and output csv file from arg
    profiler_file = sys.argv[1]
    csv_file = profiler_file + '.csv'

    # load and process JSON
    df = load_profiler_data(profiler_file, 1)
    
    # output to CSV
    csv_output(csv_file, 1, df)


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
  python3 ./dg_asm_bench_postprocess.py <profiler_output.json>
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

