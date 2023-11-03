"""
Script selects assembly data from profiler output and converts its to csv table.
Usage:

    python3 profiler_to_csv.py <profiler_file>
"""

import sys
import os
import json
import csv
from pathlib import Path
import numpy as np
import pandas as pd
import zipfile



class ProcessTag:
    """
    Collect tag processing functions.
    """

    @staticmethod
    def assembly(json_node, ph, df):
        """
        Process one assembly (Mass, Stiffness or Sources) and all subtags 
        """

        asm_name = json_node['tag']
        sum_subtags_time = sum(ch['cumul-time-sum'] for ch in json_node['children'])

        asm_time = json_node['cumul-time-sum']
        reminder_fraction = (asm_time - sum_subtags_time) / asm_time
        asm_row_data = ph.form_row(
            assembly_class=asm_name,
            time_fraction_of_reminder=reminder_fraction
        )
        return ph.df_append(df, asm_row_data)


    @staticmethod
    def asm_child(json_node, ph, df):
        integrals = {
            "assemble_volume_integrals": "cell",
            "assemble_fluxes_boundary":  "boundary_side",
            "assemble_fluxes_elem_elem": "edge",
            "assemble_fluxes_elem_side": "dimjoin"
        }

        tag = json_node['tag']
        integral = integrals.get(tag, '')
        row_data = ph.form_row(
            assembly_class=ph.parent_node['tag'],
            integral_type=integral,
        )
        return ph.df_append(df, row_data)
            

    @staticmethod
    def empty(json_node, ph, df):
        """
        We can skip some tags which are not to be processed (e.g. 'ZERO-TIME STEP')
        """
        return df


    @staticmethod
    def full_mesh(json_node, ph, df):
        """
        Process 'full_mesh' tag
        """
        mesh_repeats = json_node['call-count-sum']
        n_mesh_elements = json_node['children'][1]['call-count-sum']

        shape, dim, uniform = ph.node_path[-4]['tag'].split('_')
        size = ph.node_path[-3]['tag']
        asm_type, field_type = ph.node_path[-2]['tag'].split('_')
        ph.set_mesh_dict(
            branch=ph.program_branch, 
            commit=ph.program_revision,
            run_id=ph.run_id,
            domain_shape=shape,
            mesh_size=size,
            spacedim=dim,
            uniformity=uniform,
            n_mesh_elements=(n_mesh_elements/mesh_repeats),
            n_repeats=mesh_repeats,
            field_variant=field_type,
            assembly_variant=asm_type,
            assembly_class='',
            tag='',
            number_of_calls=0,
            integral_type='',
            time=0.0,
            time_fraction=0.0,
            time_fraction_of_reminder=np.NaN,
            parent_tag='',
            level=0
            )
        mesh_row_data = ph.form_row()
        return ph.df_append(df, mesh_row_data)


"""
Hold processed profiler data and helper variables. 
Object is used during processing of profiler JSON input.
"""
class ProfilerHandler:
    def __init__(self, program_branch, program_revision, run_id):

        # Define name of columns 
        self.column_names = ['branch', 'commit', 'run_id', 'domain_shape', 'mesh_size', 'spacedim', 'uniformity', 'n_mesh_elements', 'n_repeats', 'field_variant',   
            'assembly_variant',
                             'assembly_class', 'tag', 'number_of_calls', 'integral_type', 'time', 'time_fraction', 'time_fraction_of_reminder']
        
        self.assembly_tags = ['MassAssembly', 'StiffnessAssembly', 'SourcesAssembly']
        
        # Dictionary holds pairs tag-method, these methods are called for selected tags during processing of profiler tree
        self.tag_method = {
            'full_mesh':                     ProcessTag.full_mesh, 
            'BaseMeshReader - mesh factory': ProcessTag.empty, 
            'ZERO-TIME STEP':                ProcessTag.empty}
        self.tag_method.update({t: ProcessTag.assembly for t in self.assembly_tags})
        
        # Create other data members using during profiler processing
        self.node_path = []
        self.mesh_dict = None
        self.program_branch = program_branch
        self.program_revision = program_revision
        self.run_id = run_id


    def set_mesh_dict(self, **kwargs):
        self.mesh_dict = kwargs

    def form_row(self, **kwargs):
        time = self.current_node['cumul-time-sum']

        row = self.mesh_dict.copy()
        row.update(dict(
            tag=self.current_node['tag'],
            number_of_calls=self.current_node['call-count-sum'],
            time=time,
            time_fraction=time / self.parent_node['cumul-time-sum'],
            parent_tag=self.parent_node['tag'],
            level=len(self.node_path)
        ))
        row.update(kwargs)
        return row

    def df_append(self, df, row_df:dict):
        #row_df = pd.DataFrame(row_data, index=[0])

        if df is None:
            df = {k:[v] for k, v in row_df.items()}
            return df
        # For unknown reason df.append is deprecated.
#        c_row_df = list(row_df.columns)
#        c_df = list(df.columns)
#        assert c_row_df == c_df, f"\nrow_df: {c_row_df}\ndf:{c_df}"

        for k,v in row_df.items():
            df[k].append(v)
        return df
    @property
    def current_node(self):        
        return self.node_path[-1]
    
    @property
    def parent_node(self):
        try:
            return self.node_path[-2]
        except IndexError:
            return None
    
    def process_node_item(self, df):
        """
        Detect nodes for which we generate row  in df.
        """
        try:
            process_method = self.tag_method[self.current_node['tag']] 
        except KeyError:
            if self.parent_node is not None and self.parent_node['tag'] in self.assembly_tags:
                process_method = ProcessTag.asm_child
            else:
                return df
        return process_method(self.current_node, self, df)
    
    
    # Set base program parameters


    def location(self):
        """
        Form current node address in terms of tags.
        :return: str
        """
        return ' > '.join([node['tag'] for node in self.node_path])

    def child_nodes(self):
        if self.current_node['tag'] == "ZERO-TIME STEP":
            return []
        return self.current_node.get('children', [])
def process_node(json_node, ph, df):
    """
    Process one tag and call recursivelly processing of children
    """
    ph.node_path.append(json_node)
    #print(ph.location())
    df = ph.process_node_item(df)

    for i in ph.child_nodes():
        df = process_node(i, ph, df)
    ph.node_path.pop()
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

def profiler_path(profiler_files):
    prof_path = Path(profiler_files)
    if prof_path.suffix == '.zip':
        # path in Zip archive
        path = zipfile.Path(zipfile.ZipFile(prof_path)).iterdir()
        basename = prof_path.stem
    else:
        # path from glob pattern
        path = prof_path.parent.glob(prof_path.name)
        basename = prof_path.stem.rstrip("_*")
    return path, basename

"""
Load data from set of profiler JSON files.

Name of files must be in format '<profiler_file>_<n>.json' where:
    <profiler_file> Same for all files in set, it is given by arg 'profiler_file'
    <n> Order of file in set 1,2,3,... Number of files is given by arg 'n_runs'
"""
def load_profiler_data_(file_gen):
    # Create ProfilerHandler and DataFrame
    df = None # directory - list representation of the data table
    for run_id, file in enumerate(file_gen):
        with file.open() as f_in:
            profiler_data = json.load(f_in)
        
        # Only one descendant 'Whole Program' of root tree exists 
        whole_program = profiler_data['children'][0]
        program_branch = profiler_data['program-branch']
        program_revision = profiler_data['program-revision']
        ph = ProfilerHandler(program_branch, program_revision, run_id)
    
        df = process_node(whole_program, ph, df)
        run_id += 1
    df = pd.DataFrame(df)
    df = unify_df_values(df)
    return df

def load_profiler_data(profiler_pattern_or_zip):
    file_gen, basename = profiler_path(profiler_pattern_or_zip)
    return load_profiler_data_(file_gen)

# Perform outout to CSV file
def csv_output(csv_file, df):
    df.to_csv(path_or_buf=csv_file, header=True, index=False, mode='w')


def main():
    pass
    
    # store names of input zip archive and output csv file from arg
    profiler_pattern_or_zip = sys.argv[1]
    # make faile generator and basename
    file_gen, basename = profiler_path(profiler_pattern_or_zip)
    df =  load_profiler_data_(file_gen)

    # output to CSV
    csv_file = basename + '.csv'
    csv_output(csv_file, df)


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

