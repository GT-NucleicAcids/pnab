import json

import numpy as np


def get_results():
    """
    Read energy_data.csv file and convert to dictionary format.
    """

    data = np.genfromtxt('energy_data.csv', delimiter=',', skip_header=1)
    if data.ndim == 1:
        data = data.reshape((1, len(data)))

    keys = [i.strip() for i in open('energy_data.csv').readline().strip('\n').split(',')][1:]

    results = {}

    for conformer in range(data.shape[0]):
        results['conformer_%i.pdb' %data[conformer, 0]] = {}
        for i in range(1, data.shape[1]):
            results['conformer_%i.pdb' %data[conformer, 0]][keys[i-1]] = data[conformer, i]

    return results


def make_json(results, output_file_name='output.json'):
    """
    Outputs results in json format
    """

    with open(output_file_name, 'w') as output_file:
        json.dumps(results, output_file, indent=2)
