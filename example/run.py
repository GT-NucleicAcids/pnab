import bind

hp = bind.HelicalParameters()
hp.tilt = hp.roll = hp.shift = hp.slide = hp.buckle = hp.propeller = hp.opening = hp.shear = hp.stretch = hp.stagger = [0]
hp.twist = [32.73]
hp.rise = [2.81]
hp.inclination = [15.76]
hp.tip = [7.35]
hp.x_displacement = [-4.61]
hp.y_displacement = [-0.19]

rp = bind.RuntimeParameters()
rp.num_steps = 360000
rp.type = 'GAFF'
rp.algorithm = 'WMC'
rp.energy_filter = [10000.0, 1000000000000000.0, 1000000000000000.0, 1000000000000000.0, 1000000000000000.0]
rp.max_distance = 0.15
rp.strand = ['Adenine', 'Adenine', 'Adenine']
rp.is_double_stranded = False
rp.base_to_backbone_bond_length = 0.15

base = bind.Base()
base.file_path = 'adenine.pdb'
base.code='ADA'
base.linker = [5, 11]
base.name = 'Adenine'
base.pair_name = 'Uracil'
bases = [base]

backbone = bind.Backbone()
backbone.file_path = 'rna_bb.pdb'
backbone.interconnects = [9, 14]
backbone.linker = [12, 13]

bind.run(rp, backbone, bases, hp)
