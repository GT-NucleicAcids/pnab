from __future__ import division, absolute_import, print_function


def test_binder_1():
    """
    test binder

    Test exposed C++ classes attributes and their types
    """
    from pnab import bind

    backbone = bind.Backbone()
    backbone_attr = {'file_path': str, 'interconnects': list, 'linker': list}
    assert all([i in backbone.__dir__() for i in backbone_attr])
    assert all([type(backbone.__getattribute__(k)) is val
                for k, val in backbone_attr.items()])

    base = bind.Base()
    base_attr = {'name': str, 'code': str, 'pair_name': str, 'file_path': str, 'linker': list}
    assert all([i in base.__dir__() for i in base_attr]) 
    assert all([type(base.__getattribute__(k)) is val
                for k, val in base_attr.items()])

    runtime_parameters = bind.RuntimeParameters()
    runtime_parameters_attr = {'energy_filter': list, 'max_distance': float, 'ff_type': str, 'seed': int,
                               'search_algorithm': str, 'num_steps': int, 'dihedral_step': float,
                               'strand': list, 'is_double_stranded': bool} 
    assert all([i in runtime_parameters.__dir__() for i in runtime_parameters_attr])
    assert all([type(runtime_parameters.__getattribute__(k)) is val
                for k, val in runtime_parameters_attr.items()])

    helical_parameters = bind.HelicalParameters()
    helical_parameters_attr = ['h_twist', 'tip', 'inclination',
                               'h_rise', 'x_displacement', 'y_displacement']

    assert all([i in helical_parameters.__dir__() for i in helical_parameters_attr])
    assert all([type(helical_parameters.__getattribute__(i)) is float
                for i in helical_parameters_attr])


def test_binder_2():
    """
    test binder 2

    Test that the C++ binder runs
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
    rp.search_algorithm = 'weighted random search'
    rp.num_steps = 1000000
    rp.ff_type = 'GAFF'
    rp.energy_filter = [10000.0, 10000.0, 10000.0, 10000.0]
    rp.max_distance = 0.05
    rp.strand = ['Adenine']*5
    rp.is_double_stranded = False

    output = bind.run(rp, backbone, bases, hp, 'test')
    print(output)
