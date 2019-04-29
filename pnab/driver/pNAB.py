"""proto-Nucleic Acid Builder
proto-Nucleic Acid Builder

Main file
"""

from __future__ import division, absolute_import, print_function

import yaml
import itertools
from io import StringIO
import multiprocessing as mp
import warnings
warnings.filterwarnings("ignore", category=UserWarning)
import datetime
import numpy as np

from pnab import bind
from pnab.driver import options
try:
    from pnab.driver import widgets
except ImportError:
    pass
from pnab.driver import draw


class pNAB(object):
    """ pNAN run class

    A class that contains methods to create a pNAB run. It makes input,
    runs it, and parses results.
    """

    def __init__(self, file_path=None):
        """The constructor"""

        if file_path is not None:
            self.options = yaml.load(open(file_path, 'r'), yaml.FullLoader)
            options._validate_all_options(self.options)
            self._input_file = file_path

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
                self.options = yaml.load(open(file_path, 'r'), yaml.FullLoader)
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

        if result.size == 0:
            result = None
        elif result.ndim == 1:
            result = result.reshape(1, len(result))

        header = ''.join(['%s=%.2f, ' %(k, val) for k, val in zip(self.options['HelicalParameters'], config)])
        header = header.strip(', ')

        return [prefix, header, result] 


    def _single_result(self, results):
        """ Extract results for each helical configuration as it finishe."""

        with open('prefix.dat', 'ab') as f:
            f.write(str.encode(yaml.dump({results[0]:results[1]})))

        if results[2] is None:
            return

        header = results[1] + '\n'
        header += ('Prefix, Conformer Index, Energy (kcal/mol), Distance (A), Bond Energy, Angle Energy,' +
                  ' Torsion Energy, VDW Energy, Total Torsion Energy, RMSD (A)')

        with open('results.dat', 'ab') as f:
            np.savetxt(f, results[2], fmt='%12i'*2 + '%12.4f'*8, header=header)


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
                f.write(yaml.dump(self.options))

        config = itertools.product(*[np.random.uniform(val[0], val[1], val[2])
                                       for val in self.options['HelicalParameters'].values()])
        num_config = np.prod([val[2] for val in self.options['HelicalParameters'].values()])
        prefix = (str(i) for i in range(1, num_config + 1))

        pool = mp.Pool(mp.cpu_count(), maxtasksperchild=1)

        time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        with open('results.dat', 'w') as f: f.write('# ' + time + '\n')
        with open('prefix.dat', 'w') as f: f.write('# ' + time + '\n')

        for results in pool.imap(self._run, zip(config, prefix)):
            self._single_result(results)

        pool.close()

        if self._input_file is None:
            self._options = temp


    def get_results(self):
        """Extract the results from the run and report it to the user."""

        header = ('Prefix, Conformer Index, Energy (kcal/mol), Distance (A), Bond Energy, Angle Energy,' +
                  ' Torsion Energy, VDW Energy, Total Torsion Energy, RMSD (A)')

        time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        self.results = np.loadtxt('results.dat')
        if self.results.size == 0:
            print("No candidate found")
            return

        elif self.results.ndim == 1:
            self.results = self.results.reshape(1, len(self.results))

        self.prefix = yaml.load(open('prefix.dat'), yaml.FullLoader)

        summary = self.results[self.results[:, 2].argsort()][:10]

        np.savetxt('summary.dat', summary, fmt='%12i'*2 + '%12.4f'*8,
                   header=time + '\n' + header)

        print("There are %i candidates:\n" %len(self.results))
        print('Showing the best %i candidates ...\n' %len(summary))

        for conformer in summary:
            print(int(conformer[1]))
            print(self.prefix['%i' %conformer[0]])

            for i in range(2, len(conformer) - 1):
                print(header.split(', ')[i] + ': ' + str(conformer[i]))
            print('\n')

            draw.view_py3dmol(str(int(conformer[0])) + '_' + str(int(conformer[1])) + '.pdb')
