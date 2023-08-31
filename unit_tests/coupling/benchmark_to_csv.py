"""
Script selects assembly data from profiler output and converts its to csv table.
Usage:

    python benchmark_to_csv.py <profiler_file> <csv_file> <output_type> <mesh_size> <field_variant>
                                                           w/a           S/M/L       model/constant
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
    

# Load JSON file
profiler_file = '../' + sys.argv[1]
f_in = open(profiler_file)
profiler_data = json.load(f_in)
f_in.close()

# Only one descendant 'Whole Program' of root tree exists 
whole_program = profiler_data['children'][0]
whole_program_time = whole_program['cumul-time-sum']
print("time:", whole_program_time)

# Extract zero time step and simulation step
for i in whole_program['children']:
    if i['tag'] == "ZERO-TIME STEP": zero_step = i
    if i['tag'] == "SIMULATION-ONE STEP": simulation_step = i

zero_step_data = step_profiler_data(zero_step['children'])
simulation_data = step_profiler_data(simulation_step['children'])

print(zero_step_data)
print("------------")
print(simulation_data)


# Print header to CSV file
def write_header(writer):
    header = ['domain_shape', 'uniformity', 'spacedim', 'mesh_size', 'field_variant', 'assembly_variant', 
              'assembly_class', 'tag', 'integral_type', 'time', 'time_fraction', 'time_fraction_of_reminder']
    writer.writerow(header)

"""
Return test parameters in list (data of first columns is same in all rows)
TODO should be create dynamicaly by type of mesh and test
     domain shape Lshape, square
     uniformity: uniform, refined
     space dimension: 2, 3
     mesh size: S ~3000el, M~30k el, L~ 300k el
     field variant: model, constant
     assembly variant: empty, eval, assemble
"""
def test_params(mesh_size, field_var):
    list = ['square', 'uniform', '2', mesh_size, field_var, 'assemble']
    return list

# Print assembly data to CSV file
def asm_output_in(asm_name, assembly_data, writer, mesh_size, field_var, whole_program_time):
    integrals = {
        "assemble_volume_integrals": "cell_integral",
        "assemble_fluxes_boundary":  "boundary_side_integral",
        "assemble_fluxes_elem_elem": "edge_integral",
        "assemble_fluxes_elem_side": "dimjoin_intergral"
    }
    # Time of whole assembly
    asm_time = assembly_data['asm_time']
    asm_row = test_params(mesh_size, field_var)
    asm_row.extend([asm_name, asm_name, None, asm_time, 0, (asm_time/whole_program_time)])
    writer.writerow(asm_row)
    # Times of subprocesses
    sub_names = assembly_data.keys()
    for name in sub_names:
        if name=='asm_time': continue
        sub_time = assembly_data[name]
        sub_row = test_params(mesh_size, field_var)
        integral = None
        if name in integrals:
            integral = integrals[name]
        sub_row.extend([asm_name, name, integral, sub_time, (sub_time/asm_time), 0])
        writer.writerow(sub_row)

# Prepare print of assembly data to CSV file
def asm_output(assembly_data, writer, mesh_size, field_var, whole_program_time):
    asm_names = assembly_data.keys()
    for name in asm_names:
        asm_output_in(name, assembly_data[name], writer, mesh_size, field_var, whole_program_time)


# perform output to csv file
csv_file = '../' + sys.argv[2]
out_mode = sys.argv[3]
mesh_size = sys.argv[4]
field_var = sys.argv[5]
f_out = open(csv_file, out_mode)
writer = csv.writer(f_out)

# print header if file output mode is rewrite
if out_mode == 'w': write_header(writer)

# print whole_program data
whole_program_row = test_params(mesh_size, field_var)
whole_program_row.extend([None, whole_program['tag'], None, whole_program_time, 0, 0])
writer.writerow(whole_program_row)

# print data of assembly calls
asm_output(zero_step_data, writer, mesh_size, field_var, whole_program_time)

# close csv
f_out.close()

