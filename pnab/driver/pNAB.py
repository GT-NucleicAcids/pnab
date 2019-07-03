"""proto-Nucleic Acid Builder
proto-Nucleic Acid Builder

Main file
"""

from __future__ import division, absolute_import, print_function

import os
import yaml
import itertools
from io import StringIO
import multiprocessing as mp
import warnings
warnings.filterwarnings("ignore", category=UserWarning)
import datetime
import numpy as np

from pnab import __file__ as pnab_dir
from pnab import bind
from pnab.driver import options
try:
    from pnab.driver import jupyter_widgets
except ImportError:
    print("Running the program in the command line")

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
                  "1. Specify your options below\n"
                  "2. Use an existing input file\n")

            input_method = input('Enter input method (1, 2) [1]: ') or '1'
            input_method = input_method.strip()

            if input_method not in ['1', '2']:
                raise Exception('Please choose either "1" or "2"')

            if input_method == '1':
                self._options = options.set_options()
                self.options = None
                self._input_file = None

            else:
                file_path = input('Enter path to input [options.yaml]: ') or 'options.yaml'
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

        py_bases = [self.options[i] for i in self.options if 'Base' in i]
        bases = [bind.Base() for i in range(len(py_bases))]
        for i, b in enumerate(py_bases):
            [bases[i].__setattr__(k, val) for k, val in b.items()]

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

        with open('prefix.yaml', 'ab') as f:
            f.write(str.encode(yaml.dump({results[0]:results[1]})))

        if results[2] is None:
            return

        header = results[1] + '\n'
        header += ('Prefix, Conformer Index, Energy (kcal/mol), Distance (A), Bond Energy, Angle Energy,' +
                  ' Torsion Energy, VDW Energy, Total Torsion Energy, RMSD (A)')

        with open('results.csv', 'ab') as f:
            np.savetxt(f, results[2], delimiter=',', header=header)


    def run(self):
        """ Fuction to prepare helical configurations and run them in parallel."""

        if self._input_file is None:
            # Used Jupyter notebook; extract options
            self.options = jupyter_widgets.extract_options(self._options)
            options._validate_all_options(self.options)
            # Workaround as widget states are not serializable and cannot be used with multiprocess
            temp = self._options
            del self._options

            with open('options.yaml', 'w') as f:
                f.write(yaml.dump(self.options))

        # Add library of bases
        data_dir = os.path.join(os.path.dirname(pnab_dir), 'data')
        bases_lib = yaml.load(open(os.path.join(data_dir, 'bases_library.yaml'), 'r'), yaml.FullLoader)
        for b in bases_lib.values():
            b['file_path'] = os.path.join(data_dir, b['file_path'])
        if self.options['RuntimeParameters'].pop('pair_A_U', None):
            bases_lib['Base A']['pair_name'] = 'U'
        self.options.update(bases_lib)

        # Extract configurations
        config = itertools.product(*[np.random.uniform(val[0], val[1], val[2])
                                       for val in self.options['HelicalParameters'].values()])
        num_config = np.prod([val[2] for val in self.options['HelicalParameters'].values()])
        prefix = (str(i) for i in range(1, num_config + 1))

        # Catch interruption
        import signal
        def init_worker():
            signal.signal(signal.SIGINT, signal.SIG_IGN)

        pool = mp.Pool(mp.cpu_count(), init_worker)

        time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        with open('results.csv', 'w') as f: f.write('# ' + time + '\n')
        with open('prefix.yaml', 'w') as f: f.write('# ' + time + '\n')

        try:
            for results in pool.imap(self._run, zip(config, prefix)):
                self._single_result(results)
        except KeyboardInterrupt:
            print("Caught interruption; stopping ...")
            pool.terminate()

        pool.close()

        if self._input_file is None:
            self._options = temp


    def get_results(self):
        """Extract the results from the run and report it to the user."""

        self.prefix = yaml.load(open('prefix.yaml'), yaml.FullLoader)

        header = ('Prefix, Conformer Index, Energy (kcal/mol), Distance (A), Bond Energy, Angle Energy,' +
                  ' Torsion Energy, VDW Energy, Total Torsion Energy, RMSD (A)')

        time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        self.results = np.loadtxt('results.csv', delimiter=',')
        if self.results.size == 0:
            print("No candidate found")
            return

        elif self.results.ndim == 1:
            self.results = self.results.reshape(1, len(self.results))

        summary = self.results[self.results[:, 2].argsort()][:10]

        np.savetxt('summary.csv', summary, delimiter=',', header=time + '\n' + header)

        print("There are %i candidates:\n" %len(self.results))
        print('Showing the best %i candidates ...\n' %len(summary))

        for conformer in summary:
            print("Prefix: %i" %conformer[0])
            print(int(conformer[1]))
            print(self.prefix['%i' %conformer[0]])

            for i in range(2, len(conformer) - 1):
                print(header.split(', ')[i] + ': ' + str(conformer[i]))
            print('\n')

            draw.view_py3dmol(str(int(conformer[0])) + '_' + str(int(conformer[1])) + '.pdb')
