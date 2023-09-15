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
def load_data(profiler_file):
    # Load JSON file
    with open(profiler_file) as f_in:
        profiler_data = json.load(f_in)

    # Only one descendant 'Whole Program' of root tree exists 
    whole_program = profiler_data['children'][0]

    # Store data of all meshes to one dictionary
    prof_out_data = dict()
    # Sum time of one mesh is store to separate dictionary
    mesh_full_data = dict()

    # Iterate over data of meshes
    for i in whole_program['children']:
        key_mesh = i['tag']
        mesh_time = i['cumul-time-sum']
        mesh_repeats = i['children'][0]['call-count-sum']
        mesh_full_data[key_mesh] = (mesh_time, mesh_repeats)

        # Extract type of simulation and type of field
        for j in i['children']:
            key_type = j['tag']
            # Extract zero time step and simulation step
            for k in j['children']:
                if k['tag'] == "ZERO-TIME STEP": zero_step = k
                if k['tag'] == "SIMULATION-ONE STEP": simulation_step = k

            zero_step_data = step_profiler_data(zero_step['children'])
            simulation_data = step_profiler_data(simulation_step['children'])
            key_mesh_type = key_mesh + '_' + key_type
            prof_out_data[key_mesh_type] = simulation_data

    #print("------------")
    #print(prof_out_data)
    #print("============")
    return mesh_full_data, prof_out_data


# Print header to CSV file
def write_header(writer):
    header = ['commit', 'run_id', 'domain_shape', 'uniformity', 'spacedim', 'mesh_size', 'n_repeats', 'field_variant',  
              'assembly_variant', 'assembly_class', 'tag', 'integral_type', 'time', 'time_fraction', 'time_fraction_of_reminder']
    writer.writerow(header)


# Print data of full meshes to CSV file
def write_full_meshes(program_params, mesh_full_data, writer):
    key_names = mesh_full_data.keys()
    for key in key_names:
        mesh_split = key.split("_")
        row = [program_params[0], program_params[1], mesh_split[0], mesh_split[3], mesh_split[1], mesh_split[2], mesh_full_data[key][1], '',  
               '', '', 'full_mesh', '', mesh_full_data[key][0], '', '']
        writer.writerow(row)
    

"""
Return test parameters in list (data of first columns is same in all rows of one mesh)
Columns specifying mesh are following:
     domain shape: cube, Lshape, square
     uniformity: uniform, refined
     space dimension: 2D, 3D
     mesh size: small ~3000el, medium ~30k el, big ~300k el
     n_repeats: 10 in case of small mesh, 1 in other cases
     field variant: model, const
     assembly variant: FullAssembly, ComputeLocal, EvalField
"""
def test_params(program_params, mesh_name, mesh_full_data):
    # posible format: square_2D_small_uniform_FullAssembly_const
    separator = '_'
    params = mesh_name.split(separator)
    mesh_short_name = separator.join(params[0:4])
    list = [ program_params[0], program_params[1], params[0], params[3], params[1], params[2], mesh_full_data[mesh_short_name][1], params[5], params[4]]
    return list


# Print assembly data to CSV file
def asm_output_in(program_params, mesh_name, asm_name, assembly_data, mesh_full_data, writer):
    integrals = {
        "assemble_volume_integrals": "cell_integral",
        "assemble_fluxes_boundary":  "boundary_side_integral",
        "assemble_fluxes_elem_elem": "edge_integral",
        "assemble_fluxes_elem_side": "dimjoin_intergral"
    }
    # Time of whole assembly
    asm_time = assembly_data['asm_time']
    sum_sub_times = 0.0
    for name in assembly_data.keys():
        if name=='asm_time': continue
        sum_sub_times += assembly_data[name]
    asm_row = test_params(program_params, mesh_name, mesh_full_data)
    asm_row.extend([asm_name, asm_name, None, asm_time, 0, ((asm_time-sum_sub_times)/asm_time)])
    writer.writerow(asm_row)
    # Times of subprocesses
    for name in assembly_data.keys():
        if name=='asm_time': continue
        sub_time = assembly_data[name]
        sub_row = test_params(program_params, mesh_name, mesh_full_data)
        integral = None
        if name in integrals:
            integral = integrals[name]
        sub_row.extend([asm_name, name, integral, sub_time, (sub_time/asm_time), 0])
        writer.writerow(sub_row)


# Prepare print of mesh data to CSV file
def asm_output_mesh(program_params, mesh_name, assembly_data, mesh_full_data, writer):
    asm_names = assembly_data.keys()
    for asm_name in asm_names:
        asm_output_in(program_params, mesh_name, asm_name, assembly_data[asm_name], mesh_full_data, writer)


# Prepare print of profiler data to CSV file
def asm_output(program_params, mesh_full_data, profiler_data, writer):
    mesh_names = profiler_data.keys()
    for mesh_name in mesh_names:
        asm_output_mesh(program_params, mesh_name, profiler_data[mesh_name], mesh_full_data, writer)


# Perform outout to CSV file
def csv_output(csv_file, program_params, mesh_full_data, prof_out_data):
    out_mode = 'w'
    with open(csv_file, out_mode) as f_out:
        writer = csv.writer(f_out)

        # print header if file output mode is rewrite
        if out_mode == 'w': write_header(writer)

        # print data of full meshes
        write_full_meshes(program_params, mesh_full_data, writer)

        # print data of assembly calls
        asm_output(program_params, mesh_full_data, prof_out_data, writer)


def main():
    # store names of input json and output csv file from arg
    profiler_file = sys.argv[1]
    # commit, run_id 
    program_params = sys.argv[2:4]
    file_name_split = profiler_file.split(".")
    csv_file = file_name_split[0] + '.csv'

    # load and process JSON
    input_result = load_data(profiler_file)
    mesh_full_data = input_result[0]
    prof_out_data  = input_result[1]
    
    # output to CSV
    csv_output(csv_file, program_params, mesh_full_data, prof_out_data)



main()

