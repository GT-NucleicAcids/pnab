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

        ref_output = np.genfromtxt(os.path.join('files', f.split('.')[0] + '.csv'), delimiter=',')

        run1 = pnab.pNAB(f)
        run1.options['HelicalParameters']['is_helical'] = True
        run1.run()

        run2 = pnab.pNAB(f)
        run2.options['HelicalParameters']['is_helical'] = False
        run2.run()

        assert np.allclose(run1.results, run2.results, atol=0.4)

        if platform.system() == 'Linux' or run1.options['RuntimeParameters']['search_algorithm'] == 'systematic search':
            assert np.allclose(run1.results, ref_output, atol=0.4)
            assert np.allclose(run2.results, ref_output, atol=0.4)
