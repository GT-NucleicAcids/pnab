# coding: utf-8
import os
import sys
sys.path.insert(0, '/theoryfs2/ds/alenaizan/Gits/pNAB/driver')
import pNAB
os.chdir('/theoryfs2/ds/alenaizan/Gits/pNAB/driver/examples')


inst = pNAB.pNAB_run.init_with_questions()
print(inst.get_options())

inst.make_input('test_input.dat')
inst.make_json('test_json.json')
inst.run_input(background_run=False)
print(inst.get_results())
