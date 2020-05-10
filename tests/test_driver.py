from __future__ import division, absolute_import, print_function

import os
import platform
import glob
import numpy as np

def test_options():
    """
    test parity between the C++ attributes of classes and the python options
    """
    from pnab import bind
    from pnab.driver.options import _options_dict

    components = ['Backbone', 'Base', 'RuntimeParameters', 'HelicalParameters']
    assert all([i in _options_dict for i in components])


    assert all([i in _options_dict['Backbone'] for i in bind.Backbone.__dict__
                if not i.startswith('__')])
    assert all([i in _options_dict['Base'] for i in bind.Base.__dict__
                if not i.startswith('__')])
    assert all([i in _options_dict['RuntimeParameters'] for i in
                bind.RuntimeParameters.__dict__ if not i.startswith('__')])
    assert all([i in _options_dict['HelicalParameters'] for i in
                bind.HelicalParameters.__dict__ if not i.startswith('__')])


def test_examples():
    """
    test running with provided option file
    """
    import pnab

    examples = ['RNA.yaml', 'DNA.yaml', 'FRNA.yaml', 'LNA.yaml', 'CeNA.yaml', 'PNA.yaml', '5methylcytosine.yaml',
                'ZP.yaml', 'Hexad.yaml', 'Hexad_Antiparallel.yaml', 'adenine_cyanuric_acid.yaml']

    os.chdir(os.path.dirname(os.path.realpath(__file__)))

    for f in examples:
        print("Testing ", f)

        run = pnab.pNAB(f)
        run.run()

        if platform.system() == 'Linux' or run.options['RuntimeParameters']['search_algorithm'] == 'systematic search':
            ref_output = np.genfromtxt(os.path.join('files', f.split('.')[0] + '.csv'), delimiter=',')
            assert np.allclose(run.results, ref_output, atol=1e-4)

def test_helical_parameters():
    """
    test the equivalence between helical and step parameters.

    While the the two schemes are equivalent, numerical conversions from one scheme to the other may lead to 
    solutions that are not exactly identical. We use a loose threshold for comparison
    """

    import pnab

    examples = ['RNA.yaml', 'DNA.yaml', 'FRNA.yaml', 'LNA.yaml', 'CeNA.yaml',
                'PNA.yaml', '5methylcytosine.yaml', 'ZP.yaml']

    os.chdir(os.path.dirname(os.path.realpath(__file__)))

    for f in examples:
        print("Testing the equivalence between helical and step parameters for ", f)

        run = pnab.pNAB(f)
        run.options['RuntimeParameters']['search_algorithm'] = 'random search'
        run.options['RuntimeParameters']['num_steps'] = 1000000
        run.options['RuntimeParameters']['max_distance'] = 1
        run.options['RuntimeParameters']['energy_filter'] = [2, 4, 10, 100, 100]
        run.options['RuntimeParameters']['strand'] = ['G', 'C']
        run.options['RuntimeParameters']['num_candidates'] = 10

        run.options['HelicalParameters']['is_helical'] = True
        run.run()
        results1 = run.results

        run.options['HelicalParameters']['is_helical'] = False
        run.run()
        results2 = run.results

        assert np.allclose(results1, results2, atol=1)
