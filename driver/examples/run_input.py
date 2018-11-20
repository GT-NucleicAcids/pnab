# coding: utf-8
import os
import sys
sys.path.insert(0, '/theoryfs2/ds/alenaizan/Gits/pNAB/driver')
import pNAB
os.chdir('/theoryfs2/ds/alenaizan/Gits/pNAB/driver/examples')


inst = pNAB.pNAB_run.init_with_input('input.dat')
print(inst.get_options())

inst.make_input('test_input.dat')
inst.make_json('test_json.json')
inst.run_input(background_run=False)
print(inst.get_results())

inst2 = pNAB.pNAB_run.init_with_input('input.dat')
inst2.make_input('test_input.dat')
inst2.run_input(background_run=True)
print('see stdout.log')
print(inst2.get_run_status())
print(inst2.get_run_progress())
