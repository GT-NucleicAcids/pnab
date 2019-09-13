"""!@brief proto-Nucleic Acid Builder

This is the main file for the pnab driver. This file contains the class for calling
the pNAB C++ library. 
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

from pnab import __path__
from pnab import bind
from pnab.driver import options


class pNAB(object):
    """!@brief The proto-Nucleic Acid Builder main class

    A class that contains methods to create a pNAB run. It validates input,
    runs it, and save results.

    Examples:
    Use the code to get RNA and DNA candidates. This uses the provided example files
    in the "pnab/data" directory. The results are in three files:

    1) results.csv: Contains all the results for all the helical configurations
        requested in the run.

    2) summary.csv: Contains a summary of the best 10 results ordered by the
        total energy of the conformer.
    
    3) prefix.yaml: Containts a dictionary of the sequence of the run and the 
        corresponding helical configuration (useful when a range of values is
        given for the helical parameters).

    @code
    import pnab

    run = pnab.pNAB('RNA.yaml')
    run.run()
    @endcode

    @code
    import pnab

    run = pnab.pNAB('DNA.yaml')
    run.run()
    @endcode
    """

    def __init__(self, file_path):
        """!@brief The constructor for the pNAB run

        The constructor takes an input file as an argument and initialize a new instance.
        The input file uses a YAML format for specifying the run options. This input
        file can be generated using the graphical user interface created in @a jupyter_widgets.
        Additionally, there are a few example files in the "pnab/data" directory. This constructor
        calls @a options.validate_all_options to validate all the options for the backbone parameters,
        helical parameters, and runtime parameters. Then, it adds the library of all the defined
        nucleobases. The nucleobases are defined in the "pnab/data" directory. Whether to pair
        adenine with thymine (default) or uracil is determined here by the @a 
        pNAB.pNAB.options['RuntimeParameters']['pair_A_U'].

        The constructor creates the @a pNAB.pNAB.options attribute which contains a dictionary of all
        the defined options.

        @param file_path (str) Path to a file containing the user-defined input. This function
            tries to find the file in the current working directory. If it does not exist,
            it looks for the file in the "pnab/data" directory

        @sa jupyter_widgets
        @sa options.validate_all_options
        """

        if not os.path.isfile(file_path):
            file_path = os.path.join(__path__[0], 'data', file_path)

        ##@brief Validated options dictionary for the pNAB run
        #@sa options.validate_all_options
        self.options = yaml.load(open(file_path, 'r'), yaml.FullLoader)
        options.validate_all_options(self.options)

        # Add library of bases
        data_dir = os.path.join(__path__[0], 'data')
        bases_lib = yaml.load(open(os.path.join(data_dir, 'bases_library.yaml'), 'r'), yaml.FullLoader)
        for b in bases_lib.values():
            b['file_path'] = os.path.join(data_dir, b['file_path'])
        if self.options['RuntimeParameters'].pop('pair_A_U', None):
            bases_lib['Base A']['pair_name'] = 'U'
        self.options.update(bases_lib)

        ##@brief A string comma-separated header for the results used in the generated CSV files
        self.header = ('Prefix, Conformer Index, Distance (Angstroms), Bond Energy (kcal/mol/(nucleotide-1)), Angle Energy (kcal/mol/(nucleotide-1)), ' +
                       'Torsion Energy (kcal/mol/nucleotide), Van der Waals Energy (kcal/mol/nucleotide), ' + 
                       'Total Energy (kcal/mol/nucleotide), Fixed Torsions Energy (kcal/mol/nucleotide), ' + 
                       'RMSD relative to lowest energy conformer (Angstrom)')


    def _run(self, config):
        """!@brief Function to run one helical configuration.

        Setup the options to call the C++ code using the pybind11 bind class in @a binder.cpp.
        After the run finishes, it returns the results for processing.

        @param config (list) a list containing two entries: the list of helical parameter values for this
            configuration, and a string index for the sequnce of the run

        @returns results [prefix, header, result] A list containing the index of the run,
            a string describing the helical configuration, and a numpy array of the results

        @sa run
        @sa _single_result
        """
        config, prefix = config[0], config[1]

        # Set runtime parameters
        runtime_parameters = bind.RuntimeParameters()
        [runtime_parameters.__setattr__(k, val) for k, val in self.options['RuntimeParameters'].items()]

        # Set backbone parameters
        backbone = bind.Backbone()
        [backbone.__setattr__(k, val) for k, val in self.options['Backbone'].items()]

        # Set a list of all defined bases asd there parameters
        # We include all the bases even though not all of them are requested by the user
        py_bases = [self.options[i] for i in self.options if 'Base' in i]
        bases = [bind.Base() for i in range(len(py_bases))]
        for i, b in enumerate(py_bases):
            [bases[i].__setattr__(k, val) for k, val in b.items()]

        # Set helical parameters
        helical_parameters = bind.HelicalParameters()
        [helical_parameters.__setattr__(k, val) for k, val in zip(self.options['HelicalParameters'], config)]

        # Run code
        result = bind.run(runtime_parameters, backbone, bases, helical_parameters, prefix)

        # Get results from the comma-separated output string
        result = np.genfromtxt(StringIO(result), delimiter=',')

        if result.size == 0:
            result = None
        elif result.ndim == 1:
            result = result.reshape(1, len(result))

        # Add a header of the helical parameters
        header = ''.join(['%s=%.2f, ' %(k, val) for k, val in zip(self.options['HelicalParameters'], config)])
        header = header.strip(', ')

        results = [prefix, header, result]

        return results


    def _single_result(self, results):
        """!@brief Write results to disk for each helical configuration as it finishe.

        This is a separate function from @a pNAB.pNAB._run to allow the results to be
        written in the order that they are generated. Two files are written: "prefix.yaml"
        contains a dictionary of the sequence of the run and the corresponding helical
        configuration; "results.csv" contains the information on the accepted candidates
        for all of configurations

        @param results (list) The return value of @a pNAB.pNAB._run

        @returns None; writing to output files

        @sa _run
        @sa run
        """

        # Write prefix and helical configuration
        with open('prefix.yaml', 'ab') as f:
            f.write(str.encode(yaml.dump({results[0]:results[1]})))

        # If there are no results, return
        if results[2] is None:
            return

        # Update header
        header = results[1] + '\n'
        header += self.header

        # Write results to file
        with open('results.csv', 'ab') as f:
            np.savetxt(f, results[2], delimiter=',', header=header)


    def run(self):
        """!@brief Prepare helical configurations and run them in parallel.

        If a single value is given for a helical parameter (e.g. helical twist),
        then that value is used. If a range of values is given, then random values
        in the range will be used. The number of configurations is determined by the
        third value in @a pNAB.pNAB.options. The various helical configurations
        are run in parallel using the multiprocessing library. 

        This function writes three output files: "results.csv", "prefix.yaml",
        and "summary.csv". It renames any existing files with these names by prepending
        enough "_".

        @returns None; output files are written
        """

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

        # Rename files that have the same name
        for f in ['results.csv', 'prefix.yaml', 'summary.csv']:
            file_path = f
            while True:
                if os.path.isfile(file_path):
                    file_path = '_' + file_path
                else:
                    if os.path.isfile(f):
                        os.rename(f, file_path)
                    break

        # Write time stamps
        time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        with open('results.csv', 'w') as f: f.write('# ' + time + '\n')
        with open('prefix.yaml', 'w') as f: f.write('# ' + time + '\n')

        try:
            # Run the different helical configurations in parallel and process results
            for results in pool.imap(self._run, zip(config, prefix)):
                self._single_result(results)
        except KeyboardInterrupt:
            # If interuption is catched, terminate run and proceed
            print("Caught interruption; stopping ...")
            pool.terminate()

        pool.close()


        #Extract the results from the run and report it to the user

        ##@brief A dictionary of the sequence of the run and the 
        # corresponding helical configuration
        # For runs that have ranges of values for the helical parameters,
        # this dictionary provides a mapping between the sequence of the run
        # and the helical parameters.
        self.prefix = yaml.load(open('prefix.yaml'), yaml.FullLoader)

        # Write time stamp
        time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')

        ##@brief A numpy array containing the results of all the candidates
        #
        # The columns in the array correspond to the entries in @a pNAB.pNAB.header
        # and the helical configurations correspond to those in @a pNAB.pNAB.prefix
        self.results = np.loadtxt('results.csv', delimiter=',')

        if self.results.size == 0:
            # No results found; do not write summary file
            return

        elif self.results.ndim == 1:
            # Reshape results.
            self.results = self.results.reshape(1, len(self.results))

        # Sort by total energy
        self.results = self.results[self.results[:, 7].argsort()]

        # Save the 10 best candidates
        summary = self.results[:10]
        np.savetxt('summary.csv', summary, delimiter=',', header=time + '\n' + self.header)
