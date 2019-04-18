"""
test binder 2
"""
import os

from pNAB import bind

os.chdir(os.path.dirname(os.path.realpath(__file__)))

backbone = bind.Backbone()
backbone.file_path = 'rna_bb.pdb'
backbone.interconnects = [9, 14]
backbone.linker = [12, 13]

base = bind.Base()
base.file_path = 'adenine.pdb'
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
rp.num_steps = 1000000
rp.type = 'GAFF'
rp.algorithm = 'WMC'
rp.energy_filter = [10000.0, 10000.0, 10000.0, 10000.0, 10000.0]
rp.max_distance = 0.15
rp.strand = ['Adenine']*10
rp.is_double_stranded = False

output = bind.run(rp, backbone, bases, hp, 'test')
print(output)
