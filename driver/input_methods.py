import json
import options

"""
Methods to read and write input in input file or json formats
"""

def read_input(input_file_name):
    """
    reads input file
    """

    input_file = open(input_file_name).readlines()
    options_dict = options.options_dict()
    input_dict = {}
    num_bases = 0

    for line in input_file:
        if line.strip() and line.lstrip()[0] == '#':
            continue
        line = line.split('#')[0]

        if 'BASE PARAMETERS' in line.upper():
            num_bases += 1 
            if line.strip().upper() == 'BASE PARAMETERS':
                line = 'BASE PARAMETERS %i' %num_bases

        if line and '=' not in line:
            key1 = line.strip().upper()
            input_dict[key1] = {}

        elif '=' in line:
            input_dict[key1].update({line.split('=')[0].strip().title(): line.split('=')[1].strip()})

    if num_bases > 1:
        options.replicate_base_option(options_dict, num_bases)

    for i in options_dict:
        for j in options_dict[i]:
            try:
                options_dict[i][j].update({'new_value': input_dict[i][j]})
            except KeyError:
                options_dict[i][j].update({'new_value': options_dict[i][j]['default']})

    options.validate_options(options_dict)

    return options_dict


def read_json(json_file_name):
    """
    reads json input options.
    """

    with open(json_file_name) as f:
        json_input = json.load(f) 

    num_bases = len([i for i in json_input if 'BASE PARAMETERS' in i])
    options_dict = options.options_dict()
    if num_bases > 1:
        options.replicate_base_option(options_dict, num_bases)

    for i in options_dict:
        for j in options_dict[i]:
            try:
                options_dict[i][j].update(json_input[i][j])
            except KeyError:
                options_dict[i][j].update({'new_value': options_dict[i][j]['default']})

    options.validate_options(options_dict)
    
    return options_dict


def make_input(options_dict, input_file_name='input.dat'):
    """
    A method to generate input file from options_dict.
    Takes dictionary of options and an optional input file name.
    """
    
    input_file = open(input_file_name, 'w')

    for key1 in sorted(options_dict):
        input_file.write(key1 + '\n')

        for key2 in options_dict[key1]:
            if not isinstance(options_dict[key1][key2]['new_value'], tuple):
                input_file.write('    ' + key2 + '=' + str(options_dict[key1][key2]['new_value']) + '\n')
            else:
                input_file.write('    ' + key2 + '=')
                for elem in options_dict[key1][key2]['new_value'][:-1]:
                    input_file.write(str(elem) + ',')
                input_file.write(str(options_dict[key1][key2]['new_value'][-1]) + '\n')

        input_file.write('\n')

    input_file.close()


def make_json(options_dict, json_file_name='input.json'):
    """
    A method to generate json input file from options_dict.
    Takes dictionary of options and an optional json input file name.
    """
    
    json_input = {}
    for key1 in sorted(options_dict):
        json_input[key1] = {}

        for key2 in options_dict[key1]:
            json_input[key1][key2] = {'new_value': options_dict[key1][key2]['new_value']}

    with open(json_file_name, 'w') as input_file:
        json.dump(json_input, input_file, indent=2)

    return json_input
