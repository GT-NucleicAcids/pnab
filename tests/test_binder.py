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
                               'strand': list, 'is_hexad': bool, 'build_strand': list, 'strand_orientation': list,
                               'weighting_temperature': float, 'monte_carlo_temperature': float, 'mutation_rate': float,
                               'crossover_rate': float, 'population_size': int, 'glycosidic_bond_distance': float,
                               'num_candidates': int}
    assert all([i in runtime_parameters.__dir__() for i in runtime_parameters_attr])
    assert all([type(runtime_parameters.__getattribute__(k)) is val
                for k, val in runtime_parameters_attr.items()])

    helical_parameters = bind.HelicalParameters()
    helical_parameters_attr = ['h_twist', 'tip', 'inclination', 'h_rise', 'x_displacement', 'y_displacement',
                               'twist', 'roll', 'tilt', 'rise', 'slide', 'shift',
                               'buckle', 'propeller', 'opening', 'stretch', 'shear', 'stagger']

    assert all([i in helical_parameters.__dir__() for i in helical_parameters_attr])
    assert all([type(helical_parameters.__getattribute__(i)) is float
                for i in helical_parameters_attr])
    assert type(helical_parameters.__getattribute__('is_helical')) is bool


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
    base.code = 'A'
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
    rp.energy_filter = [2.0, 5.0, 10.0, 10000.0, 10000.0]
    rp.max_distance = 0.2
    rp.strand = ['Adenine']*5
    rp.num_candidates = 2

    output = bind.run(rp, backbone, bases, hp, 'test')
    print(output)

def test_helical_parameters():
    """
    Test the equivalence between helical and step parameters.

    The test compares 3 sets of values for helical and step parameters.
    Using the two schemes should give identical results for the accepted candidates.
    To ensure correct comparisons, the step parameter values are set to many digits.
    """
    import os

    from pnab import bind

    os.chdir(os.path.dirname(os.path.realpath(__file__)))

    parameters = [
                  [[1, 2, 3, 4, 5, 6], [-0.4691332223, 0.3125034304, 2.9558478770, -0.5222738076, 0.4178190461, 5.9626387714]],
                  [[-1, -2, -3, -4, -5, -6], [-0.4691332223, 0.3125034304, -2.9558478770, -0.5222738076, 0.4178190461, -5.9626387714]],
                  [[-2, 5, 3, -15, 20, -40], [2.2177307312, 0.4662275467, 4.1753711904, 13.2973947701, 9.9730460776, -36.5122480218]]
                 ]

    backbone = bind.Backbone()
    backbone.file_path = os.path.join('files', 'rna_bb.pdb')
    backbone.interconnects = [10, 1]
    backbone.linker = [13, 14]

    base = bind.Base()
    base.file_path = os.path.join('files', 'adenine.pdb')
    base.code = 'A'
    base.linker = [5, 11]
    base.name = 'Adenine'
    base.pair_name = 'Uracil'
    bases = [base]

    rp = bind.RuntimeParameters()
    rp.search_algorithm = 'weighted random search'
    rp.num_steps = 10
    rp.ff_type = 'GAFF'
    rp.energy_filter = [1e10]*5
    rp.max_distance = 1e10
    rp.strand = ['Adenine']*5
    rp.num_candidates = 10

    for p in parameters:
        hp = bind.HelicalParameters()
        hp.is_helical = True
        hp.x_displacement = p[0][0]
        hp.y_displacement = p[0][1]
        hp.h_rise = p[0][2]
        hp.inclination = p[0][3]
        hp.tip = p[0][4]
        hp.h_twist = p[0][5]

        output_helical = bind.run(rp, backbone, bases, hp, 'test')

        hp.is_helical = False
        hp.shift = p[1][0]
        hp.slide = p[1][1]
        hp.rise = p[1][2]
        hp.tilt = p[1][3]
        hp.roll = p[1][4]
        hp.twist = p[1][5]

        output_step = bind.run(rp, backbone, bases, hp, 'test')

        assert output_helical == output_step
