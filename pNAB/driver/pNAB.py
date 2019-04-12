"""proto-Nucleic Acid Builder
proto-Nucleic Acid Builder

details
"""

from __future__ import division, absolute_import, print_function

import json
import itertools
from io import StringIO
import multiprocessing as mp
import warnings
warnings.filterwarnings("ignore", category=UserWarning)
import datetime
import numpy as np

from pNAB import bind
from pNAB.driver import options
from pNAB.driver import widgets
from pNAB.driver import draw


class pNAB(object):
    """ pNAN run class

    A class that contains methods to create a pNAB run. It makes input,
    runs it, and parses results.
    """

    def __init__(self, file_path=None):
        """The constructor"""

        if file_path is not None:
            self.options = json.loads(open(file_path, 'r').read())
            options._validate_all_options(self.options)

        else:
            print("There are two methods to use the Nucleic Acid Builder:\n"
                  "1. Answer a list of questions\n"
                  "2. Use an existing input file\n")

            input_method = input('Enter input method (1, 2) [1]: ') or '1'

            if input_method not in ['1', '2']:
                raise Exception('Please choose either "1" or "2"')

            if input_method == '1':
                self._options = options.set_options()
                self.options = None
                self._input_file = None

            else:
                file_path = input('Enter path to input [options.dat]: ') or 'options.dat'
                self.options = json.loads(open(file_path, 'r').read())
                options._validate_all_options(self.options)
                self._input_file = file_path


    def _run(self, config):
        """ function to run one helical configuration."""
        config, prefix = config[0], config[1]

        runtime_parameters = bind.RuntimeParameters()
        [runtime_parameters.__setattr__(k, val) for k, val in self.options['RuntimeParameters'].items()]

        backbone = bind.Backbone()
        [backbone.__setattr__(k, val) for k, val in self.options['Backbone'].items()]

        num_bases = len([i for i in self.options if 'Base' in i])
        bases = [bind.Base() for i in range(num_bases)]
        for i in range(len(bases)):
            [bases[i].__setattr__(k, val) for k, val in self.options['Base %i' %(i + 1)].items()]


        helical_parameters = bind.HelicalParameters()
        [helical_parameters.__setattr__(k, val) for k, val in zip(self.options['HelicalParameters'], config)]
        result = bind.run(runtime_parameters, backbone, bases, helical_parameters, prefix)
        result = np.genfromtxt(StringIO(result), delimiter=',')

        if result.ndim == 1:
            result = result.reshape(1, len(result))

        if result.shape[1] != 0:
            header = ''.join(['%s=%.2f, ' %(k, val) for k, val in zip(self.options['HelicalParameters'], config)])
            header = header.strip(', ')
        else:
            header = result = None

        return [prefix, header, result] 


    def run(self):
        """ Fuction to prepare helical configurations and run them in parallel."""

        if self._input_file is None:
            # Used Jupyter notebook; extract options
            self.options = widgets.extract_options(self._options)
            options._validate_all_options(self.options)
            # Workaround as widget states are not serializable and cannot be used with multiprocess
            temp = self._options
            del self._options
            with open('options.dat', 'w') as f:
                f.write(json.dumps(self.options, indent=4))

        config = list(itertools.product(*[np.random.uniform(val[0], val[1], val[2]) 
                                       for k, val in self.options['HelicalParameters'].items()]))

        prefix = [str(i) for i in range(1, len(config) + 1)]

        pool = mp.Pool(mp.cpu_count())

        self._results = pool.map(self._run, zip(config, prefix))

        pool.close()

        if self._input_file is None:
            self._options = temp


    def get_results(self):
        """Extract the results from the run and report it to the user."""

        header = ('Conformer Index, Energy (kcal/mol), Distance (A), Bond Energy, Angle Energy,' +
                  ' Torsion Energy, VDW Energy, Total Torsion Energy, RMSD (A)')

        results = self._results
        f = open('results.dat', 'ab')
        config_dict = {}
        prefix_dict = {}
        time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        all_results = np.zeros((0, 10))
        for i in range(len(results)):
            if results[i][1] is None:
                continue
            prefix_dict[str(int(results[i][0]))] = results[i][1]
            prefix_str = 'prefix: ' + results[i][0] + '\n'
            np.savetxt(f, results[i][2], fmt='%12i' + '%12.2f'*8, header=time + '\n' + prefix_str + results[i][1] + '\n' + header)
            for conf in results[i][2][:, 0]:
                config_dict[str(int(results[i][0])) + '_' + str(int(conf))] = (prefix_str + results[i][1]).strip(',')
            prefix = np.array([int(results[i][0])]*results[i][2].shape[0]).reshape((results[i][2].shape[0], 1))
            all_results = np.vstack((all_results, np.hstack((prefix, results[i][2]))))
        f.close()

        all_results = all_results[all_results[:, 2].argsort()]
        self.config_dict = config_dict
        self.results = all_results

        max_show = 10 if len(self.results) > 10 else len(self.results)
        f = open('summary.dat', 'ab')
        np.savetxt(f, self.results[:max_show], fmt='%12i' + '%12i' + '%12.2f'*8,
                   header=time + '\n'  + 'Prefix, ' + header )
        f.close()

        f = open('prefix.dat', 'ab')
        f.write(str.encode(time + '\n'))
        for k, v in prefix_dict.items():
            f.write(str.encode(k + ': ' + v + '\n'))
        f.close()

        print("There are %i candidates:\n" %len(self.results))
        print('Showing the best %i candidates ...\n' %max_show)

        for conformer in list(self.results)[:max_show]:
            print('\n\n' + str(int(conformer[1])))
            print('\n' +  config_dict[str(int(conformer[0])) + '_' + str(int(conformer[1]))]) 

            for i in range(2, len(conformer)):
                print(header.split(', ')[i-1] + ': ' + str(conformer[i]))

            draw.view_py3dmol(str(int(conformer[0])) + '_' + str(int(conformer[1])) + '.pdb')
