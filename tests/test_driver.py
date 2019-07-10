from __future__ import division, absolute_import, print_function


def test_options():
    """
    test parity between the C++ attributes of classes and the python options
    """
    from pnab import bind
    from pnab.driver.options import _options_dict

    components = ['Backbone', 'Base 1', 'RuntimeParameters', 'HelicalParameters']
    assert all([i in _options_dict for i in components])


    assert all([i in _options_dict['Backbone'] for i in bind.Backbone.__dict__
                if not i.startswith('__')])
    assert all([i in _options_dict['Base 1'] for i in bind.Base.__dict__
                if not i.startswith('__')])
    assert all([i in _options_dict['RuntimeParameters'] for i in
                bind.RuntimeParameters.__dict__ if not i.startswith('__')])
    assert all([i in _options_dict['HelicalParameters'] for i in
                bind.HelicalParameters.__dict__ if not i.startswith('__')])


def test_run():
    """
    test running with provided option file
    """
    import os

    import pnab

    os.chdir(os.path.dirname(os.path.realpath(__file__)))

    run = pnab.pNAB('files/options_rna.yaml')
    run.run()
    run.get_results()


def test_run_range():
    """
    test running with provided option file
    """
    import os

    import pnab

    os.chdir(os.path.dirname(os.path.realpath(__file__)))

    run = pnab.pNAB('files/options_rna_range.yaml')
    run.run()
    run.get_results()

    # Confirm the number of tested configurations
    assert len(run.prefix) == 15


def test_duplex():
    """
    test generating duplex with provided option file
    """
    import os

    import pnab

    os.chdir(os.path.dirname(os.path.realpath(__file__)))

    run = pnab.pNAB('files/options_dna.yaml')
    run.run()
    run.get_results()


def test_hexad():
    """
    test generating hexad with provided option file
    """
    import os

    import pnab

    os.chdir(os.path.dirname(os.path.realpath(__file__)))

    run = pnab.pNAB('files/options_hexad.yaml')
    run.run()
    run.get_results()
