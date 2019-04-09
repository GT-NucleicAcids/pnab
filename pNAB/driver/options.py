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

from pNAB.driver import draw

def set_options():
    """A method to get user-defined options."""

    options_dict = _ask_options_questions()

    return options_dict


def _ask_options_questions():
    """A method that ask the user for options."""

    options_dict = {}
    num_bases = int(input("Enter number of distinct bases [1]: ") or 1)
    _replicate_base_option(num_bases)

    for k1 in sorted(_options_dict):
        print('\n', k1)
        options_dict[k1] = {}

        for k2 in _options_dict[k1]:
            options_dict[k1][k2] = (input("    Enter %s [%s]: " %(_options_dict[k1][k2]['glossory'], 
                                                                  _options_dict[k1][k2]['default']))
                                                               or _options_dict[k1][k2]['default'])

            options_dict[k1][k2] = _options_dict[k1][k2]['validation'](
                                                        options_dict[k1][k2])

            if k2 == 'file_path':
                draw.view_py3dmol(options_dict[k1][k2], label=True)

            if str(options_dict[k1][k2]) == 'None':
                raise Exception('Value of %s in %s has not been set correctly' %(k2, k1))

    return options_dict


def _validate_all_options(options):
    """A method to validate all options.

    Used when the user provides an input file
    """
    for k1 in options:
        k1_val = k1
        if 'Base' in k1:
            k1_val = 'Base 1'

        for k2 in options[k1]:
            options[k1][k2] = _options_dict[k1_val][k2]['validation'](options[k1][k2])

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
                                                             'glossory': ('path to the file containing the molecular' +
                                                                         ' structure of the backbone'),
                                                             'default': 'backbone.pdb',
                                                             'validation': lambda x: _validate_input_file(x),
                                                             }

_options_dict['Backbone']['interconnects'] = {
                                                        'glossory': 'two atom indices that connect to the two backbones',
                                                        'default': None,
                                                        'validation': lambda x: _validate_atom_indices(x),
                                                        }

_options_dict['Backbone']['linker'] = {
                                                       'glossory': ('two atom indices that form the vector' + 
                                                                   ' connecting backbone to base molecule'),
                                                       'default': None,
                                                       'validation': lambda x: _validate_atom_indices(x),
                                                       }

# Base Parameters
_options_dict['Base 1'] = {}

_options_dict['Base 1']['file_path'] = {
                                                       'glossory': 'path to file containing the structure of the base',
                                                       'default': 'base.pdb',
                                                       'validation': lambda x: _validate_input_file(x),
                                                       }
_options_dict['Base 1']['linker'] = {
                                                         'glossory': ('two atoms that define a vector connecting to backbone' +
                                                                     ' atoms that connect to the base'),
                                                         'default': None,
                                                         'validation': lambda x: _validate_atom_indices(x),
                                                         }

_options_dict['Base 1']['code'] = {
                                             'glossory': 'three-letter code used to identify residue in PDB format',
                                             'default': 'RES',
                                             'validation': lambda x: str(x),
                                             }

_options_dict['Base 1']['name'] = {
                                             'glossory': 'base name',
                                             'default': 'base',
                                             'validation': lambda x: str(x),
                                             }

_options_dict['Base 1']['pair_name'] = {
                                                  'glossory': 'name of the pairing base',
                                                  'default': 'pairing_base',
                                                  'validation': lambda x: str(x),
                                                  }

# Helical Parameters
_options_dict['HelicalParameters'] = {}
#_options_dict['HelicalParameters']['tilt'] = {
#                                              'glossory': 'tilt',
#                                              'default': 0.0,
#                                              'validation': lambda x: _validate_helical_parameters(x),
#                                              }
#
#_options_dict['HelicalParameters']['roll'] = {
#                                              'glossory': 'roll',
#                                              'default': 0.0,
#                                              'validation': lambda x: _validate_helical_parameters(x),
#                                              }

_options_dict['HelicalParameters']['twist'] = {
                                               'glossory': 'twist',
                                               'default': 0.0,
                                               'validation': lambda x: _validate_helical_parameters(x),
                                               }

_options_dict['HelicalParameters']['shift'] = {
                                               'glossory': 'shift',
                                               'default': 0.0,
                                               'validation': lambda x: _validate_helical_parameters(x),
                                               }

_options_dict['HelicalParameters']['slide'] = {
                                               'glossory': 'slide',
                                               'default': 0.0,
                                               'validation': lambda x: _validate_helical_parameters(x),
                                               }

_options_dict['HelicalParameters']['rise'] = {
                                              'glossory': 'rise',
                                              'default': 0.0,
                                              'validation': lambda x: _validate_helical_parameters(x),
                                              }

#_options_dict['HelicalParameters']['buckle'] = {
#                                                'glossory': 'buckle',
#                                                'default': 0.0,
#                                                'validation': lambda x: _validate_helical_parameters(x),
#                                                }
#
#_options_dict['HelicalParameters']['propeller'] = {
#                                                   'glossory': 'propeller',
#                                                   'default': 0.0,
#                                                   'validation': lambda x: _validate_helical_parameters(x),
#                                                   }
#
#_options_dict['HelicalParameters']['opening'] = {
#                                                 'glossory': 'opening',
#                                                 'default': 0.0,
#                                                 'validation': lambda x: _validate_helical_parameters(x),
#                                                 }
#
#_options_dict['HelicalParameters']['shear'] = {
#                                               'glossory': 'shear',
#                                               'default': 0.0,
#                                               'validation': lambda x: _validate_helical_parameters(x),
#                                               }
#
#_options_dict['HelicalParameters']['stretch'] = {
#                                                 'glossory': 'stretch',
#                                                 'default': 0.0,
#                                                 'validation': lambda x: _validate_helical_parameters(x),
#                                                 }
#
#_options_dict['HelicalParameters']['stagger'] = {
#                                                 'glossory': 'stagger',
#                                                 'default': 0.0,
#                                                 'validation': lambda x: _validate_helical_parameters(x),
#                                                 }

_options_dict['HelicalParameters']['inclination'] = {
                                                     'glossory': 'inclination',
                                                     'default': 0.0,
                                                     'validation': lambda x: _validate_helical_parameters(x),
                                                     }

_options_dict['HelicalParameters']['tip'] = {
                                             'glossory': 'tip',
                                             'default': 0.0,
                                             'validation': lambda x: _validate_helical_parameters(x),
                                             }

_options_dict['HelicalParameters']['x_displacement'] = {
                                                        'glossory': 'x-displacement',
                                                        'default': 0.0,
                                                        'validation': lambda x: _validate_helical_parameters(x),
                                                        }

_options_dict['HelicalParameters']['y_displacement'] = {
                                                        'glossory': 'y-displacement',
                                                        'default': 0.0,
                                                        'validation': lambda x: _validate_helical_parameters(x),
                                                        }

# Runtime Parameters
_options_dict['RuntimeParameters'] = {}
_options_dict['RuntimeParameters']['num_steps'] = {
                                                   'glossory': 'number of points to search',
                                                   'default': 1000000,
                                                   'validation': lambda x: int(x),
                                                   }

_options_dict['RuntimeParameters']['type'] = {
                                              'glossory': 'force field type',
                                              'default': 'GAFF',
                                              'validation': lambda x: str(x).upper(),
                                              }

# Not used yet.
#_options_dict['RuntimeParameters']['parameter_file'] = {
#                                                        'glossory': 'additional parameter file',
#                                                        'default': '',
#                                                        'validation': lambda x: _validate_input_file(x),
#                                                        }

_options_dict['RuntimeParameters']['algorithm'] = {
                                                   'glossory': 'three letter code for search algorithm',
                                                   'default': 'WMC',
                                                   'validation': lambda x: _validate_algorithm(x),
                                                   }

_options_dict['RuntimeParameters']['energy_filter'] = {
                                                       'glossory': ('list of energy filters (kcal/mol): \n\n' +
                                                                    '1. cuttoff for the total energy of the entire molecule' + 
                                                                    ' resulting from sequence specified\n' + 
                                                                    '2. cuttoff for energy from angles between three consecutive atoms\n' +
                                                                    '3. cuttoff for energy from bond between newly bonded atoms\n' +
                                                                    '4. cuttoff for total van der Waals energy\n' + 
                                                                    '5. cuttoff for energy from the new torsion formed' + 
                                                                    ' bewteen linking backbones\n'),
                                                       'default': (1e5, 1e5, 1e5, 1e5, 1e5),
                                                       'validation': lambda x: _validate_energy_filter(x), 
                                                      }

_options_dict['RuntimeParameters']['max_distance'] = {
                                                      'glossory': ('the maximum distance between atom' +
                                                                   ' linkers in backbone. Conformers that' + 
                                                                   ' do not meet this criteria are' + 
                                                                   ' automatically rejected and do not even' + 
                                                                   ' pass to energy analysis where a chain is built'),
                                                      'default': 0.5,
                                                      'validation': lambda x: float(x),
                                                      }

_options_dict['RuntimeParameters']['strand'] = {
                                                'glossory': 'list of base names in strand',
                                                'default': None,
                                                'validation': lambda x: ([i.strip() for i in x.split(',') if i] if isinstance(x, str)
                                                                         else [str(i).strip() for i in tuple(x) if i]),
                                                }

_options_dict['RuntimeParameters']['is_double_stranded'] = {
                                                            'glossory': 'whether to generate double stranded conformers',
                                                            'default': False,
                                                            'validation': lambda x: bool(eval(x.title())) if isinstance(x, str) else bool(x),
                                                            }
