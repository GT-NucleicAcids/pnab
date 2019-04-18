"""options
file for preparing and validation options

This is an option file that defines a dictionary containing all the available options for the builder code.
options_dict has at least four keys: Backbone, Base, HelicalParameters, and RuntimeParameters.
Note that there might be more than one Base keys that are numbered to account for
the presence of multiple bases. Each base has its own associated options.
Each of these contains a dictionary of available options in each category. Each dictionary contains a glossory
that describe the option, a default value, a validation scheme, a new value, and validation
confirmation.
"""

from __future__ import division, absolute_import, print_function

import os
import copy

try:
    from pNAB.driver import widgets
except:
    pass

def set_options():
    """A method to get user-defined options."""

    options_dict = _ask_options_questions()

    return options_dict


def _ask_options_questions():
    """A method that ask the user for options."""

    options_dict = {}
    options_dict = widgets.display_widgets(_options_dict)

    return options_dict


def _validate_all_options(options):
    """A method to validate all options.

    Used when the user provides an input file
    """

    num_bases = len(['Base' for i in options if 'Base' in i])
    _replicate_base_option(num_bases)

    for k1 in _options_dict:
        for k2 in _options_dict[k1]:
            options[k1][k2] = _options_dict[k1][k2]['validation'](options[k1][k2])


def _replicate_base_option(num_bases):
    """replicates initial options for bases if the number of distinct bases is greater than 1."""

    for i in range(2, num_bases + 1):
        _options_dict.update({'Base %i' %i: copy.deepcopy(_options_dict['Base 1'])})


def _validate_input_file(file_name):
    """Method to validate provided path to a geomerty file. Check whether the file exists or not."""

    file_name = str(file_name)
    if not os.path.isfile(file_name):
        raise Exception("Cannot find file: %s" %file_name)

    return file_name


def _validate_atom_indices(x):
    """Method to validate provided tuple of atom indices."""

    x = list(eval(x)) if isinstance(x, str) else list(x)
    if len(x) != 2:
        raise Exception("Incorrect number of atoms. Must provide two indices.")

    for i in x:
        if i == 0:
            raise Exception("Use 1-based index")

    return x


def _validate_helical_parameters(x):
    """Method to validate provided helical parameters."""

    if isinstance(x, str):
        x = eval(x)

    if isinstance(x, (float, int)):
        x = [x, x, 1]
    else:
        x = list(x)

    if len(x) == 1:
        x = [x[0], x[0], 1]
    elif len(x) == 2:
        # Default of 5 values tested for a range
        x.append(5)
    elif len(x) > 3:
        raise Exception("Helical parameters must have at most three values")

    return x


def _validate_energy_filter(x):
    """Method to validate provided energy filter."""

    if isinstance(x, str):
        x = eval(x)

    x = list(x)
    if len(x) != 5:
        raise Exception("Five values must be provided for the energy filter")

    for i in x:
        if not isinstance(i, (float, int)):
            raise Exception("Provide valid numbers")

    return x
   
def _validate_algorithm(algorithm):
    """method to validate algorithm. Checks if the requested algorith is available."""

    available_algorithms = ['WMC']
    algorithm = str(algorithm).upper()
    if algorithm in available_algorithms:
        return algorithm
    else:
        raise Exception('Available alogrithms are %s' %str(available_algorithms))


# Set glossory of options, default values and validation methods

_options_dict = {}

# Backbone Parameter
_options_dict['Backbone'] = {}
_options_dict['Backbone']['file_path'] = {
                                                             'glossory': ('Path to the file containing the molecular' +
                                                                         ' structure of the backbone'),
                                                             'default': 'backbone.pdb',
                                                             'validation': lambda x: _validate_input_file(x),
                                                             }

_options_dict['Backbone']['interconnects'] = {
                                                        'glossory': 'Two atoms connecting to the two backbones',
                                                        'default': None,
                                                        'validation': lambda x: _validate_atom_indices(x),
                                                        }

_options_dict['Backbone']['linker'] = {
                                                       'glossory': ('Two atoms forming the vector' + 
                                                                   ' connecting to base'),
                                                       'default': None,
                                                       'validation': lambda x: _validate_atom_indices(x),
                                                       }

# Base Parameters
_options_dict['Base 1'] = {}

_options_dict['Base 1']['file_path'] = {
                                                       'glossory': 'Path to file containing the structure of the base',
                                                       'default': 'base.pdb',
                                                       'validation': lambda x: _validate_input_file(x),
                                                       }
_options_dict['Base 1']['linker'] = {
                                                         'glossory': ('Two atoms forming a vector connecting to backbone'),
                                                         'default': None,
                                                         'validation': lambda x: _validate_atom_indices(x),
                                                         }

_options_dict['Base 1']['code'] = {
                                             'glossory': 'Three-letter code',
                                             'default': 'RES',
                                             'validation': lambda x: str(x),
                                             }

_options_dict['Base 1']['name'] = {
                                             'glossory': 'Base name',
                                             'default': 'base',
                                             'validation': lambda x: str(x),
                                             }

_options_dict['Base 1']['pair_name'] = {
                                                  'glossory': 'Name of the pairing base',
                                                  'default': 'pairing_base',
                                                  'validation': lambda x: str(x),
                                                  }

# Helical Parameters
_options_dict['HelicalParameters'] = {}
#_options_dict['HelicalParameters']['tilt'] = {
#                                              'glossory': 'Tilt',
#                                              'default': 0.0,
#                                              'validation': lambda x: _validate_helical_parameters(x),
#                                              }
#
#_options_dict['HelicalParameters']['roll'] = {
#                                              'glossory': 'Roll',
#                                              'default': 0.0,
#                                              'validation': lambda x: _validate_helical_parameters(x),
#                                              }

_options_dict['HelicalParameters']['twist'] = {
                                               'glossory': 'Twist',
                                               'default': 0.0,
                                               'validation': lambda x: _validate_helical_parameters(x),
                                               }

_options_dict['HelicalParameters']['shift'] = {
                                               'glossory': 'Shift',
                                               'default': 0.0,
                                               'validation': lambda x: _validate_helical_parameters(x),
                                               }

_options_dict['HelicalParameters']['slide'] = {
                                               'glossory': 'Slide',
                                               'default': 0.0,
                                               'validation': lambda x: _validate_helical_parameters(x),
                                               }

_options_dict['HelicalParameters']['rise'] = {
                                              'glossory': 'Rise',
                                              'default': 0.0,
                                              'validation': lambda x: _validate_helical_parameters(x),
                                              }

#_options_dict['HelicalParameters']['buckle'] = {
#                                                'glossory': 'Buckle',
#                                                'default': 0.0,
#                                                'validation': lambda x: _validate_helical_parameters(x),
#                                                }
#
#_options_dict['HelicalParameters']['propeller'] = {
#                                                   'glossory': 'Propeller',
#                                                   'default': 0.0,
#                                                   'validation': lambda x: _validate_helical_parameters(x),
#                                                   }
#
#_options_dict['HelicalParameters']['opening'] = {
#                                                 'glossory': 'Opening',
#                                                 'default': 0.0,
#                                                 'validation': lambda x: _validate_helical_parameters(x),
#                                                 }
#
#_options_dict['HelicalParameters']['shear'] = {
#                                               'glossory': 'Shear',
#                                               'default': 0.0,
#                                               'validation': lambda x: _validate_helical_parameters(x),
#                                               }
#
#_options_dict['HelicalParameters']['stretch'] = {
#                                                 'glossory': 'Stretch',
#                                                 'default': 0.0,
#                                                 'validation': lambda x: _validate_helical_parameters(x),
#                                                 }
#
#_options_dict['HelicalParameters']['stagger'] = {
#                                                 'glossory': 'Stagger',
#                                                 'default': 0.0,
#                                                 'validation': lambda x: _validate_helical_parameters(x),
#                                                 }

_options_dict['HelicalParameters']['inclination'] = {
                                                     'glossory': 'Inclination',
                                                     'default': 0.0,
                                                     'validation': lambda x: _validate_helical_parameters(x),
                                                     }

_options_dict['HelicalParameters']['tip'] = {
                                             'glossory': 'Tip',
                                             'default': 0.0,
                                             'validation': lambda x: _validate_helical_parameters(x),
                                             }

_options_dict['HelicalParameters']['x_displacement'] = {
                                                        'glossory': 'X-Displacement',
                                                        'default': 0.0,
                                                        'validation': lambda x: _validate_helical_parameters(x),
                                                        }

_options_dict['HelicalParameters']['y_displacement'] = {
                                                        'glossory': 'Y-Displacement',
                                                        'default': 0.0,
                                                        'validation': lambda x: _validate_helical_parameters(x),
                                                        }

# Runtime Parameters
_options_dict['RuntimeParameters'] = {}
_options_dict['RuntimeParameters']['num_steps'] = {
                                                   'glossory': 'Number of points to search',
                                                   'default': 1000000,
                                                   'validation': lambda x: int(x),
                                                   }

_options_dict['RuntimeParameters']['type'] = {
                                              'glossory': 'Force field type',
                                              'default': 'GAFF',
                                              'validation': lambda x: str(x).upper(),
                                              }

# Not used yet.
#_options_dict['RuntimeParameters']['parameter_file'] = {
#                                                        'glossory': 'Additional parameter file',
#                                                        'default': '',
#                                                        'validation': lambda x: _validate_input_file(x),
#                                                        }

_options_dict['RuntimeParameters']['algorithm'] = {
                                                   'glossory': 'Three letter code for search algorithm',
                                                   'default': 'WMC',
                                                   'validation': lambda x: _validate_algorithm(x),
                                                   }

_options_dict['RuntimeParameters']['energy_filter'] = {
                                                       'glossory': ('Cuttoff for the total energy of the molecule\n' + 
                                                                    'Cuttoff for energy from angles between three consecutive atoms\n' +
                                                                    'Cuttoff for energy from bond between newly bonded atoms\n' +
                                                                    'Cuttoff for total van der Waals energy\n' + 
                                                                    'Cuttoff for torsional energy between backbone linkers\n'),
                                                       'default': (1e10, 5, 5, 0, 1e0),
                                                       'validation': lambda x: _validate_energy_filter(x), 
                                                      }

_options_dict['RuntimeParameters']['max_distance'] = {
                                                      'glossory': ('The maximum distance between atom' +
                                                                   ' linkers in backbone.' 
                                                                   ),
                                                      'default': 0.05,
                                                      'validation': lambda x: float(x),
                                                      }

_options_dict['RuntimeParameters']['strand'] = {
                                                'glossory': 'List of base names in strand',
                                                'default': None,
                                                'validation': lambda x: ([i.strip() for i in x.split(',') if i] if isinstance(x, str)
                                                                         else [str(i).strip() for i in tuple(x) if i]),
                                                }

_options_dict['RuntimeParameters']['is_double_stranded'] = {
                                                            'glossory': 'Double strands',
                                                            'default': False,
                                                            'validation': lambda x: bool(eval(x.title())) if isinstance(x, str) else bool(x),
                                                            }
