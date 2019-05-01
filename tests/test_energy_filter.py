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
    backbone.interconnects = [9, 14]
    backbone.linker = [12, 13]

    base = bind.Base()
    base.file_path = 'files/adenine.pdb'
    base.code = 'ADA'
    base.linker = [5, 11]
    base.name = 'Adenine'
    base.pair_name = 'Uracil'
    bases = [base]

    hp = bind.HelicalParameters()
    hp.shift = 0
    hp.slide = 0
    hp.twist = 32.73
    hp.rise = 2.81
    hp.inclination = 15.76
    hp.tip = 7.35
    hp.x_displacement = -4.61
    hp.y_displacement = -0.19

    rp = bind.RuntimeParameters()
    rp.num_steps = 1
    rp.type = 'GAFF'
    rp.algorithm = 'WMC'
    rp.energy_filter = [1e100]*5
    rp.max_distance = 1e100
    rp.strand = ['Adenine']*5
    rp.is_double_stranded = False

    output1 = bind.run(rp, backbone, bases, hp, '1')
    output1 = np.genfromtxt(StringIO(output1), delimiter=',')

    rp.num_steps = 100000000
    rp.energy_filter = [output1[2], output1[5], output1[4], output1[7], output1[8]]
    rp.max_distance = 0.2

    output2 = bind.run(rp, backbone, bases, hp, '2')
    output2 = np.genfromtxt(StringIO(output2), delimiter=',')

    if output2.size == 0:
        return
    elif output2.ndim == 1:
        output2 = output2.reshape(1, len(output2))

    filter_index = [2, 4, 5, 7, 8]
    for i in filter_index:
        assert all([val < output1[i] for val in output2[:, i]])
