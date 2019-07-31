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

from pnab import __path__


def _validate_all_options(options):
    """A method to validate all options.

    Used when the user provides an input file
    """

    for k1 in _options_dict:
        if k1 not in options:
            options['k1'] = {}
        for k2 in _options_dict[k1]:
            if k2 not in options[k1]:
                options[k1][k2] = _options_dict[k1][k2]['default']
            options[k1][k2] = _options_dict[k1][k2]['validation'](options[k1][k2])

    if options['RuntimeParameters']['is_double_stranded'] and ("X" in options['RuntimeParameters']['strand'] or "Y" in options['RuntimeParameters']['strand']):
        raise Exception("Cannot build double strands for triaminopyrimidine or cyanuric acid")

    elif options['RuntimeParameters']['is_hexad']:
        for i in ['A', 'G', 'C', 'U', 'T']:
            if i in options['RuntimeParameters']['strand']:
                raise Exception("Cannot build hexads for canonical nucleobases")


def _validate_input_file(file_name):
    """Method to validate provided path to a geomerty file. Check whether the file exists or not."""

    file_name = str(file_name)
    if not os.path.isfile(file_name):
        file_name = os.path.join(__path__[0], 'data', file_name)
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
    if len(x) != 4:
        raise Exception("Four values must be provided for the energy filter")

    for i in x:
        if not isinstance(i, (float, int)):
            raise Exception("Provide valid numbers")

    return x
   
def _validate_base_name(name):
    name = str(name)
    if len(name) > 1:
        raise Exception("Base name must be one letter") 

    return name

def _validate_strand(strand):
    strand = list(strand)
    if "X" in strand or "Y" in strand:
        for i in ["A", "G", "C", "T", "U"]:
            if i in strand:
                raise Exception("Cannot combine canonincal and non-canonical nucleobases")

    return strand

def _validate_strand_orientation(strand):
    strand = list(strand)
    for i, val in enumerate(strand):
        strand[i] = bool(eval(val.title())) if isinstance(val, str) else bool(val)

    return strand




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
                                               'default': [0.0, 0.0, 1],
                                               'validation': lambda x: _validate_helical_parameters(x),
                                               }

_options_dict['HelicalParameters']['inclination'] = {
                                                     'glossory': 'Inclination',
                                                     'default': [0.0, 0.0, 1],
                                                     'validation': lambda x: _validate_helical_parameters(x),
                                                     }

_options_dict['HelicalParameters']['tip'] = {
                                             'glossory': 'Tip',
                                             'default': [0.0, 0.0, 1],
                                             'validation': lambda x: _validate_helical_parameters(x),
                                             }

_options_dict['HelicalParameters']['rise'] = {
                                              'glossory': 'Rise',
                                              'default': [0.0, 0.0, 1],
                                              'validation': lambda x: _validate_helical_parameters(x),
                                              }

_options_dict['HelicalParameters']['x_displacement'] = {
                                                        'glossory': 'X-Displacement',
                                                        'default': [0.0, 0.0, 1],
                                                        'validation': lambda x: _validate_helical_parameters(x),
                                                        }

_options_dict['HelicalParameters']['y_displacement'] = {
                                                        'glossory': 'Y-Displacement',
                                                        'default': [0.0, 0.0, 1],
                                                        'validation': lambda x: _validate_helical_parameters(x),
                                                        }


_options_dict['HelicalParameters']['shift'] = {
                                               'glossory': 'Shift',
                                               'default': [0.0, 0.0, 1],
                                               'validation': lambda x: _validate_helical_parameters(x),
                                               }

_options_dict['HelicalParameters']['slide'] = {
                                               'glossory': 'Slide',
                                               'default': [0.0, 0.0, 1],
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

# Runtime Parameters
_options_dict['RuntimeParameters'] = {}
_options_dict['RuntimeParameters']['num_steps'] = {
                                                   'glossory': 'Number of points to search over dihedral angles',
                                                   'default': 100000000,
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

_options_dict['RuntimeParameters']['energy_filter'] = {
                                                       'glossory': ('Maximum energy per bond for newly formed bonds in the backbone\n' +
                                                                    'Maximum energy per agnle for newly formed angles in the backbone\n' +
                                                                    'Maximum van der Waals energy per nucleotide\n' +
                                                                    'Maximum total energy per nucleotide\n'),
                                                       'default': (3, 3, 0, 10000000000),
                                                       'validation': lambda x: _validate_energy_filter(x), 
                                                      }

_options_dict['RuntimeParameters']['max_distance'] = {
                                                      'glossory': ('The maximum distance between atom' +
                                                                   ' linkers in backbone' 
                                                                   ),
                                                      'default': 0.05,
                                                      'validation': lambda x: float(x),
                                                      }

_options_dict['RuntimeParameters']['strand'] = {
                                                'glossory': 'FASTA string for nucleotide sequence (e.g. GCAT or XYXY) ',
                                                'default': None,
                                                'validation': lambda x: _validate_strand(x),
                                                }

_options_dict['RuntimeParameters']['is_double_stranded'] = {
                                                            'glossory': 'Double strands',
                                                            'default': False,
                                                            'validation': lambda x: bool(eval(x.title())) if isinstance(x, str) else bool(x),
                                                            }
_options_dict['RuntimeParameters']['pair_A_U'] = {
                                                  'glossory': 'Pair A with U (Default is A-T pairing)',
                                                  'default': False,
                                                  'validation': lambda x: bool(eval(x.title())) if isinstance(x, str) else bool(x),
                                                  }
_options_dict['RuntimeParameters']['is_hexad'] = {
                                                  'glossory': 'Hexad strands',
                                                  'default': False,
                                                  'validation': lambda x: bool(eval(x.title())) if isinstance(x, str) else bool(x),
                                                  }
_options_dict['RuntimeParameters']['strand_orientation'] = {
                                                  'glossory': 'Orientation of each strand in the hexad (up or down)',
                                                  'default': [True, True, True, True, True, True],
                                                  'validation': lambda x: _validate_strand_orientation(x),
                                                  }

