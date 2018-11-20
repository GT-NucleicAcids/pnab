import subprocess
import glob
import copy

import numpy as np

import options
import input_methods
import results


class pNAB_run(object):
    """
    A class that contains methods to create a pNAB run. It parses input,
    runs input, and parses results.
    """

    default_name = 'default'

    def __init__(self, name=default_name):
        self.name = name
        self._options = None
        self._run_process = None
        self._run_successfully = None
        self._input_file_name = None

    @classmethod
    def init_with_questions(cls, name=default_name):
        """
        initialize a pNAB_run instance with a set of questions to the user
        """

        instance = cls()
        instance._options = options.set_options(options.options_dict())

        return instance
    

    @classmethod
    def init_with_input(cls, name=default_name, input_file_name='input.dat'):
        """
        initialize a pNAB_run instance with an input file
        """

        instance = cls()
        instance._options = input_methods.read_input(input_file_name)
        instance._input_file_name = input_file_name

        return instance


    @classmethod
    def init_with_json(cls, name=default_name, json_file_name='input.json'):
        """
        initialize a pNAB_run instance with a json input file or string
        """

        instance = cls()
        instance._options = input_methods.read_json(json_file_name)

        return instance


    def make_input(self, input_file_name=None):
        """
        makes input file
        """

        if input_file_name is None:
            input_file_name = self.name + '.dat'

        input_methods.make_input(self._options, input_file_name)
        self._input_file_name = input_file_name


    def make_json(self, json_file_path=None):
        """
        makes json input file
        """

        if json_file_path is None:
            json_file_path = self.name + '.json'

        input_methods.make_json(self._options, json_file_path)


    def run_input(self, background_run=False):
        """
        Fuction to run input and return execution status
        """

        if self._run_successfully:
            raise Exception('This pNAB instance has already been run successfully' + 
                            ' or is currently running.')

        run = ['pNAB', self._input_file_name]

        if background_run:
            with open(self.name + '.log', 'w') as stdout:
                self._run_process = subprocess.Popen(run, stdout=stdout)

        else:
            run = subprocess.run(run)
            if run.returncode != 0:
                self._run_successfully = False
                raise Exception('Program did not run successfully.')

            else:
                self._run_successfully = True


    def get_options(self):
        """
        Function to get the options associated with self
        """

        return copy.deepcopy(self._options)


    def get_run_status(self):
        """
        Function to get the run status for background jobs.
        status can be: still running, False or True
        """

        if self._run_process is not None:
            status = self._run_process.poll()
            if status is None:
                self._run_successfully = 'still running'
            elif status != 0:
                self._run_successfully = False
            else:
                self._run_successfully = True

        return self._run_successfully


    def get_run_progress(self):
        """
        Function to get the progress of the run for jobs running in the background.
        """

        if self._run_process is not None:
            stdout = open(self.name + '.log').readlines()[-3:]

            return ''.join(stdout)

        else:
            return 'See standart output.'


    def get_results(self):
        """
        Function to get results from output files
        """

        return results.get_results()
