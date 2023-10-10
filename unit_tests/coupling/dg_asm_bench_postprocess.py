"""
Script selects assembly data from profiler output and converts its to csv table.
Usage:

    python3 profiler_to_csv.py <profiler_file>
"""

import sys
import json
import csv


# Store data of one assembly to dictionary
def asm_profiler_data(json_in):
    out_data = dict()
    out_data['asm_time'] = json_in['cumul-time-sum']
    for i in json_in['children']:
        key = i['tag']
        val = i['cumul-time-sum']
        out_data[key] = val
    return out_data


# Store data of set of assembly to dictionary
def step_profiler_data(json_in):
    out_data = dict()
    for i in json_in:
        if i['tag'] == "assembly":
            for j in i['children']:
                key = j['tag']
                val = asm_profiler_data(j)
                out_data[key] = val
    return out_data


# Load data from profilerJSON file
def load_data(profiler_file, whole_program_time):
    # Load JSON file
    with open(profiler_file) as f_in:
        profiler_data = json.load(f_in)

    # Only one descendant 'Whole Program' of root tree exists 
    whole_program = profiler_data['children'][0]
    whole_program_time = whole_program['cumul-time-sum']
    #print("time:", whole_program_time)

    # Store data of all meshes to one dictionary
    prof_out_data = dict()

    # Iterate over data of meshes
    for i in whole_program['children']:
        key = i['tag']

        # Extract zero time step and simulation step
        for j in i['children']:
            if j['tag'] == "ZERO-TIME STEP": zero_step = j
            if j['tag'] == "SIMULATION-ONE STEP": simulation_step = j

        zero_step_data = step_profiler_data(zero_step['children'])
        simulation_data = step_profiler_data(simulation_step['children'])
        prof_out_data[key] = simulation_data

    #print("------------")
    #print(prof_out_data)
    #print("============")
    return prof_out_data;


# Print header to CSV file
def write_header(writer):
    header = ['domain_shape', 'uniformity', 'spacedim', 'mesh_size', 'field_variant', 'assembly_variant', 
              'assembly_class', 'tag', 'integral_type', 'time', 'time_fraction', 'time_fraction_of_reminder']
    writer.writerow(header)


"""
Return test parameters in list (data of first columns is same in all rows)
TODO should be create dynamicaly by type of mesh and test
     domain shape cube, Lshape, square
     uniformity: uniform, refined
     space dimension: 2, 3
     mesh size: S ~3000el, M~30k el, L~ 300k el
     field variant: model, constant
     assembly variant: empty, eval, assemble
"""
def test_params(mesh_name, field_var):
    # square_2D_small_uniform
    params = mesh_name.split("_")
    list = [ params[0], params[3], params[1], params[2], field_var, 'assemble']
    return list


# Print assembly data to CSV file
def asm_output_in(mesh_name, asm_name, assembly_data, writer, field_var, whole_program_time):
    integrals = {
        "assemble_volume_integrals": "cell_integral",
        "assemble_fluxes_boundary":  "boundary_side_integral",
        "assemble_fluxes_elem_elem": "edge_integral",
        "assemble_fluxes_elem_side": "dimjoin_intergral"
    }
    # Time of whole assembly
    asm_time = assembly_data['asm_time']
    asm_row = test_params(mesh_name, field_var)
    asm_row.extend([asm_name, asm_name, None, asm_time, 0, (asm_time/whole_program_time)])
    writer.writerow(asm_row)
    # Times of subprocesses
    sub_names = assembly_data.keys()
    for name in sub_names:
        if name=='asm_time': continue
        sub_time = assembly_data[name]
        sub_row = test_params(mesh_name, field_var)
        integral = None
        if name in integrals:
            integral = integrals[name]
        sub_row.extend([asm_name, name, integral, sub_time, (sub_time/asm_time), 0])
        writer.writerow(sub_row)


# Prepare print of mesh data to CSV file
def asm_output_mesh(mesh_name, assembly_data, writer, field_var, whole_program_time):
    asm_names = assembly_data.keys()
    for asm_name in asm_names:
        asm_output_in(mesh_name, asm_name, assembly_data[asm_name], writer, field_var, whole_program_time)


# Prepare print of profiler data to CSV file
def asm_output(profiler_data, writer, field_var, whole_program_time):
    mesh_names = profiler_data.keys()
    for mesh_name in mesh_names:
        asm_output_mesh(mesh_name, profiler_data[mesh_name], writer, field_var, whole_program_time)


# Perform outout to CSV file
def csv_output(csv_file, prof_out_data, whole_program_time):
    out_mode = 'w'
    field_var = 'full'
    with open(csv_file, out_mode) as f_out:
        writer = csv.writer(f_out)

        # print header if file output mode is rewrite
        if out_mode == 'w': write_header(writer)

        # print whole_program data
        # TODO whole program is not used, will be replaced by time of one mesh
        #whole_program_row = [ None, None, None, None, field_var, 'assemble']
        #whole_program_row.extend([None, whole_program['tag'], None, whole_program_time, 0, 0])
        #writer.writerow(whole_program_row)

        # print data of assembly calls
        asm_output(prof_out_data, writer, field_var, whole_program_time)


def main():
    # store names of input json and output csv file from arg
    profiler_file = sys.argv[1]
    file_name_split = profiler_file.split(".")
    csv_file = file_name_split[0] + '.csv'

    whole_program_time = 1.0         # temporary
    prof_out_data = load_data(profiler_file, whole_program_time)
    csv_output(csv_file, prof_out_data, whole_program_time)



main()

