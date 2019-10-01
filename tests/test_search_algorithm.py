from __future__ import division, absolute_import, print_function


def test_binder_2():
    """
    test binder 2

    Test that all search algorithm work
    """
    import os

    from pnab import bind

    os.chdir(os.path.dirname(os.path.realpath(__file__)))

    backbone = bind.Backbone()
    backbone.file_path = os.path.join('files', 'rna_bb.pdb')
    backbone.interconnects = [10, 1]
    backbone.linker = [13, 14]

    base = bind.Base()
    base.file_path = os.path.join('files', 'adenine.pdb')
    base.code = 'ADA'
    base.linker = [5, 11]
    base.name = 'Adenine'
    base.pair_name = 'Uracil'
    bases = [base]

    hp = bind.HelicalParameters()
    hp.h_twist = 32.39
    hp.h_rise = 2.53
    hp.inclination = 22.9
    hp.tip = 0.08
    hp.x_displacement = -4.54
    hp.y_displacement = -0.02

    rp = bind.RuntimeParameters()
    rp.num_steps = 1000000
    rp.dihedral_step = 4
    rp.ff_type = 'GAFF'
    rp.energy_filter = [10000.0, 10000.0, 10000.0, 10000.0]
    rp.max_distance = 0.05
    rp.strand = ['Adenine']*3

    for s in ['weighted monte carlo search', 'monte carlo search', 'weighted random search', 'random search', 'systematic search']:
        print(s)
        rp.search_algorithm = s
        output = bind.run(rp, backbone, bases, hp, 'test')
        print(output)
