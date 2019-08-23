from __future__ import division, absolute_import, print_function

def test_energy_filter():
    """
    test that energy filters are enforced and that the program
    generates lower energy conformers
    """
    import os
    from io import StringIO

    import numpy as np

    from pnab import bind

    os.chdir(os.path.dirname(os.path.realpath(__file__)))

    backbone = bind.Backbone()
    backbone.file_path = 'files/rna_bb.pdb'
    backbone.interconnects = [10, 1]
    backbone.linker = [13, 14]

    base = bind.Base()
    base.file_path = 'files/adenine.pdb'
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
    rp.search_algorithm = 'weighted monte carlo search'
    rp.num_steps = 1
    rp.type = 'GAFF'
    rp.energy_filter = [1e100]*5
    rp.max_distance = 1e100
    rp.strand = ['Adenine']*5
    rp.is_double_stranded = False

    output1 = bind.run(rp, backbone, bases, hp, '1')
    output1 = np.genfromtxt(StringIO(output1), delimiter=',')

    rp.num_steps = 1000000
    rp.energy_filter = [output1[3], output1[4], output1[5], output1[6], output1[7]]
    rp.max_distance = 0.05

    output2 = bind.run(rp, backbone, bases, hp, '2')
    output2 = np.genfromtxt(StringIO(output2), delimiter=',')

    if output2.size == 0:
        return
    elif output2.ndim == 1:
        output2 = output2.reshape(1, len(output2))

    filter_index = [3, 4, 5, 6, 7]
    for i in filter_index:
        assert all([val < output1[i] for val in output2[:, i]])
