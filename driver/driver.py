import time
from IPython.display import clear_output

import pNAB
import view


"""
Driver for the jupyter notebook that simplifies parsing input, running, and outputting results
"""


def initialize():
    """
    input driver function for the nucleic acid builder.
    """

    input_method = input('Enter input method (1, 2) [1]: ') or '1'

    if input_method not in ['1', '2']:
        raise Exception('Please choose either "1" or "2"')

    if input_method == '1':
        inst = pNAB.pNAB_run.init_with_questions()

    elif input_method == '2':
        input_file_path = input('Enter path to input [input.dat]: ') or 'input.dat'
        inst = pNAB.pNAB_run.init_with_input(input_file_name=input_file_path)

    return inst


def run(inst):
    """
    driver runner function
    """

    if inst._input_file_name is None:
        inst.make_input()

    inst.run_input()


def get_results(inst):
    """
    driver results function
    """

    results = inst.get_results()
    jsmol_path = input('To view the molecules, please provide path to jsmol directory [None]:') or None
    print("There are %i candidates:\n" %len(results))

    for conformer in results:
        print(conformer)

        for key in results[conformer]:
            print(key, ': ', results[conformer][key])

        if jsmol_path is not None:
            a = view.view(jsmol_path, conformer)
            a.view()
            time.sleep(15)
            clear_output()
        

    return results
