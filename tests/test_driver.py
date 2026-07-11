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
                'ZP.yaml', 'Hexad.yaml', 'Hexad_Antiparallel.yaml', 'adenine_cyanuric_acid.yaml', 'adenine_cyanuric_acid_flipped.yaml']

    if platform.system() != 'Linux': # test times in other platforms exceed allowed limit for travis CI and appveyor
        examples = examples[:7]

    os.chdir(os.path.dirname(os.path.realpath(__file__)))

    for f in examples:
        print("Testing ", f)

        run = pnab.pNAB(f)
        run.run()

        if platform.system() == 'Linux':
            ref_output = np.genfromtxt(os.path.join('files', f.split('.')[0] + '.csv'), delimiter=',')
        elif platform.system() == 'Windows':
            ref_output = np.genfromtxt(os.path.join('files', f.split('.')[0] + '_windows.csv'), delimiter=',')
        elif platform.system() == 'Darwin':
            ref_output = np.genfromtxt(os.path.join('files', f.split('.')[0] + '_mac.csv'), delimiter=',')
            
        assert np.allclose(run.results, ref_output, atol=1e-4)


def test_linkage_planarity_filter():
    """
    test the max_linkage_deviation filter for non-planar inter-residue linkages

    The peptide linkage in the PNA example sits ~43 degrees from planar. Because the
    linkage dihedral is set by the periodic build and is not sampled by the search, a
    strict max_linkage_deviation rejects every candidate, while a permissive value keeps
    the candidate and reports it in a "Linkage Dihedral" column.
    """
    import pnab

    os.chdir(os.path.dirname(os.path.realpath(__file__)))

    # Strict threshold: the non-planar linkage is filtered out, so nothing is accepted
    run = pnab.pNAB('PNA.yaml')
    run.options['RuntimeParameters']['max_linkage_deviation'] = 8.0
    run.options['RuntimeParameters']['num_steps'] = 100000
    run.run()
    assert run.results.size == 0

    # Permissive threshold: the candidate is accepted and its linkage dihedral reported
    run = pnab.pNAB('PNA.yaml')
    run.options['RuntimeParameters']['max_linkage_deviation'] = 90.0
    run.run()
    assert run.results.size != 0
    assert 'Linkage Dihedral 1 (degrees)' in run.header
