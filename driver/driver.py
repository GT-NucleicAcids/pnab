import pNAB
import draw


"""
Driver for pNAB notebook that simplifies parsing input, running, and outputting results
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
    max_show = 10 if len(results) > 10 else len(results)
    print("There are %i candidates:\n" %len(results))
    print('Showing the best %i candidates ...\n' %max_show)

    for conformer in list(results)[:max_show]:
        print('\n', conformer)

        for key in results[conformer]:
            print(key, ': ', results[conformer][key])

        draw.view_nglviewer(conformer)

    return results
