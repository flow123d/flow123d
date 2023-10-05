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
        
        # Create empty DataFrame object      
        self.df = pd.DataFrame( columns=self.column_names )
        
        # Create other data members using during profiler processing
        self.common_row_data = []
        self.profiler_subtree = []
        self.sum_subtags_time = 0.0


# Load annd store data of one assembly to DataFrame
def load_asm_data(json_in, asm_name, ph):
    integrals = {
        "assemble_volume_integrals": "cell_integral",
        "assemble_fluxes_boundary":  "boundary_side_integral",
        "assemble_fluxes_elem_elem": "edge_integral",
        "assemble_fluxes_elem_side": "dimjoin_intergral"
    }
    
    ph.profiler_subtree = json_in
    ph.sum_subtags_time = 0.0
    for i in json_in['children']:
        tag = i['tag']
        time = i['cumul-time-sum']
        integral = None
        if tag in integrals:
            integral = integrals[tag]
        tag_data = ph.common_row_data.copy()
        tag_data.extend( [asm_name, tag, i['call-count-sum'], integral, time, (time/ph.profiler_subtree['cumul-time-sum']), 0] )
        ph.df.loc[len(ph.df.index)] = tag_data
        ph.sum_subtags_time += time
    
    asm_time = json_in['cumul-time-sum']
    asm_data = ph.common_row_data.copy()
    asm_data.extend( [asm_name, asm_name, json_in['call-count-sum'], None, asm_time, 0, ((asm_time-ph.sum_subtags_time)/asm_time)] )
    ph.df.loc[len(ph.df.index)] = asm_data


# Load one simulation step, iterate over assembly sub-tags
def load_simulation_step(json_in, ph):
    for i in json_in:
        if i['tag'] == "assembly":
            for j in i['children']:
                asm_name = j['tag']
                load_asm_data(j, asm_name, ph)


# Load data from profiler JSON file
def load_profiler_data(ph):
    # Load JSON file
    with open(ph.profiler_file) as f_in:
        profiler_data = json.load(f_in)
    
    # Only one descendant 'Whole Program' of root tree exists 
    whole_program = profiler_data['children'][0]

    # Iterate over data of meshes
    for i in whole_program['children']:
        key_mesh_name = i['tag']
        for j in i['children']:
            # Process profiler data of one mesh to data table
            key_mesh_size = j['tag']
            mesh_time = j['cumul-time-sum']
            key_mesh_params = key_mesh_name.split('_')

            # Extract type of simulation and type of field
            for k in j['children']:
                key_type = k['tag']
                # Extract zero time step and simulation step
                for m in k['children']:
                    if m['tag'] == "ZERO-TIME STEP":
                        mesh_repeats = m['call-count-sum']
                    if m['tag'] == "SIMULATION-ONE STEP": 
                        simulation_step = m
                    if m['tag'] == "n_mesh_elements": 
                        n_mesh_elements = m['call-count-sum']

                type_params = key_type.split('_')
                ph.common_row_data = [ ph.program_params[0], ph.program_params[1], key_mesh_params[0], key_mesh_size, key_mesh_params[1], 
                    key_mesh_params[2], (n_mesh_elements/mesh_repeats), mesh_repeats, type_params[0], type_params[1] ]
                load_simulation_step(simulation_step['children'], ph)

            # Perform data of one mesh to data table
            ph.df.loc[len(ph.df.index)] = [ ph.program_params[0], ph.program_params[1], key_mesh_params[0], key_mesh_size, key_mesh_params[1], 
                key_mesh_params[2], (n_mesh_elements/mesh_repeats), mesh_repeats, '', '', '', 'full_mesh', j['call-count-sum']*mesh_repeats, '', 
                mesh_time, '', '' ]
            
    #print("------------")
    #print(ph.df)
    #print("============")


# Perform outout to CSV file
def csv_output(csv_file, program_params, df):
    out_mode = 'a'
    out_header = False
    if program_params[1] == '1': 
        out_mode = 'w'
        out_header = True
    df.to_csv(path_or_buf=csv_file, header=out_header, index=False, mode=out_mode)


def main():
    # store names of input json and output csv file from arg
    profiler_file = sys.argv[1]
    # commit, run_id 
    program_params = sys.argv[2:4]
    file_name_split = profiler_file.split(".")
    csv_file = file_name_split[0] + '.csv'

    # load and process JSON
    ph = ProfilerHandler(profiler_file, program_params)
    load_profiler_data(ph)
    
    # output to CSV
    csv_output(csv_file, program_params, ph.df)



main()

