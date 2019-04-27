from __future__ import division, absolute_import, print_function
import pytest


def test_options():
    """
    test parity between the C++ classes and the python options
    """
    from pnab import bind
    from pnab.driver.options import _options_dict

    components = ['Backbone', 'Base 1', 'RuntimeParameters', 'HelicalParameters']
    assert all([i in _options_dict for i in components])


    assert all([i in bind.Backbone.__dict__ for i in _options_dict['Backbone']])
    assert all([i in bind.Base.__dict__ for i in _options_dict['Base 1']])
    assert all([i in bind.RuntimeParameters.__dict__ for i in _options_dict['RuntimeParameters']])
    assert all([i in bind.HelicalParameters.__dict__ for i in _options_dict['HelicalParameters']])


def test_run():
    """
    test running with provided option file
    """
    import os

    import pnab

    os.chdir(os.path.dirname(os.path.realpath(__file__)))

    run = pnab.pNAB('files/options_rna.dat')
    run.run()
    run.get_results()
    print(run.results)


#def test_run2():
#    """
#    test running with provided option file. Run over ranges of helical parameters:
#    five twist values and three inclination values
#    """
#    import os
#
#    import pnab
#
#    os.chdir(os.path.dirname(os.path.realpath(__file__)))
#
#    run = pnab.pNAB('files/options_rna_range.dat')
#    run.run()
#    assert len(run._results) == 15
#
#    run.get_results()
#
#    for prefix, conformer in zip(run.results[:, 0], run.results[:, 1]):
#        assert '%i_%i' %(prefix, conformer) in run.config_dict
#
#    print(run.results)
