"""
This is an option file that defines a dictionary containing all the available options for the builder code.
options_dict has at least four keys: BACKBONE PARAMETERS, BASE PARAMETERS, HELICAL PARAMETERS, and RUNTIME
PARAMETERS. Note that there might be more than one BASE PARAMETERS keys that are numbered to account for
the presence of multiple bases. Each base has its own associated options.
Each of these contains a dictionary of available options in each category. Each dictionary contains a glossory
that describe the option, a default value, a validation scheme, a new value, and validation
confirmation.
"""

import os
import copy


def set_options(options_dict):
    """
    A method to get user-defined options
    """

    _ask_options_questions(options_dict)
    validate_options(options_dict)

    return options_dict


def _ask_options_questions(options_dict):
    """
    A method that ask the user for options
    """

    num_bases = int(input("Enter number of distinct bases [1]: ") or 1)
    replicate_base_option(options_dict, num_bases)

    for k1 in sorted(options_dict):
        print('\n', k1)
        for k2 in options_dict[k1]:
            options_dict[k1][k2]['new_value'] = (input("    Enter %s [%s]: " %(options_dict[k1][k2]['glossory'], 
                                                                               options_dict[k1][k2]['default']))
                                                 or options_dict[k1][k2]['default'])


def replicate_base_option(options_dict, num_bases):
    """
    replicates initial options for bases if the number of distinct bases is greater than 1
    """

    for i in range(2, num_bases + 1):
        options_dict.update({'BASE PARAMETERS %i' %i: copy.deepcopy(options_dict['BASE PARAMETERS 1'])})


def validate_options(options_dict):
    """
    validate all the options in options_dict according to the validation scheme of each option.
    Update new_value with the validated value.
    """

    for key1 in options_dict:
        for key2 in options_dict[key1]:
            options_dict[key1][key2]['new_value'] = options_dict[key1][key2]['validation'](
                                                        options_dict[key1][key2]['new_value'])
            if str(options_dict[key1][key2]['new_value']) == 'None':
                raise Exception('Value of %s in %s has not been set to a correct value' %(key2, key1))
            else:
                options_dict[key1][key2]['validated'] = True


def _validate_input_file(file_name):
    """
    Method to validate provided path to a geomerty file. Check whether the file exists or not.
    """

    file_name = str(file_name)
    os.path.isfile(file_name)
    return file_name


def _validate_atom_indicies(x):
    """
    Method to validate provided tuple of atom indices
    """

    x = tuple(eval(x)) if isinstance(x, str) else tuple(x)
    if len(x) != 2:
        raise Exception("Incorrect number of atoms. Must provide two indices.")

    return x


# Set glossory of options, default values and validation methods

_options_dict = {}

# Backbone Parameter
_options_dict['BACKBONE PARAMETERS'] = {}
_options_dict['BACKBONE PARAMETERS']['Interconnects'] = {
                                                        'glossory': 'two atom indices that connect to the two backbones',
                                                        'default': None,
                                                        'validation': lambda x: _validate_atom_indicies(x),
                                                        }

_options_dict['BACKBONE PARAMETERS']['Base_Connect'] = {
                                                       'glossory': ('two atom indicies that form the vector' + 
                                                                   ' connecting backbone to base molecule'),
                                                       'default': None,
                                                       'validation': lambda x: _validate_atom_indicies(x),
                                                       }

_options_dict['BACKBONE PARAMETERS']['Backbone_File_Path'] = {
                                                             'glossory': ('path to the file containing the molecular' +
                                                                         ' structure of the backbone'),
                                                             'default': 'backbone.pdb',
                                                             'validation': lambda x: _validate_input_file(x),
                                                             }

# Base Parameters
_options_dict['BASE PARAMETERS 1'] = {}
_options_dict['BASE PARAMETERS 1']['Backbone_Connect'] = {
                                                         'glossory': ('two atoms that define a vector connecting to backbone' +
                                                                     ' atoms that connect to the base'),
                                                         'default': None,
                                                         'validation': lambda x: _validate_atom_indicies(x),
                                                         }

_options_dict['BASE PARAMETERS 1']['Code'] = {
                                             'glossory': 'three-letter code used to identify residue in PDB format',
                                             'default': 'RES',
                                             'validation': lambda x: str(x),
                                             }

_options_dict['BASE PARAMETERS 1']['Name'] = {
                                             'glossory': 'base name',
                                             'default': 'base',
                                             'validation': lambda x: str(x),
                                             }

_options_dict['BASE PARAMETERS 1']['Base_File_Path'] = {
                                                       'glossory': 'path to file containing the structure of the base',
                                                       'default': 'base.pdb',
                                                       'validation': lambda x: _validate_input_file(x),
                                                       }

_options_dict['BASE PARAMETERS 1']['Pair_Name'] = {
                                                  'glossory': 'name of the pairing base',
                                                  'default': 'pairing_base',
                                                  'validation': lambda x: str(x),
                                                  }

# Helical Parameters
_options_dict['HELICAL PARAMETERS'] = {}
_options_dict['HELICAL PARAMETERS']['Tilt'] = {
                                              'glossory': 'tilt',
                                              'default': 0.0,
                                              'validation': lambda x: float(x),
                                              }

_options_dict['HELICAL PARAMETERS']['Roll'] = {
                                              'glossory': 'roll',
                                              'default': 0.0,
                                              'validation': lambda x: float(x),
                                              }

_options_dict['HELICAL PARAMETERS']['Twist'] = {
                                               'glossory': 'twist',
                                               'default': 0.0,
                                               'validation': lambda x: float(x),
                                               }

_options_dict['HELICAL PARAMETERS']['Shift'] = {
                                               'glossory': 'shift',
                                               'default': 0.0,
                                               'validation': lambda x: float(x),
                                               }

_options_dict['HELICAL PARAMETERS']['Slide'] = {
                                               'glossory': 'slide',
                                               'default': 0.0,
                                               'validation': lambda x: float(x),
                                               }

_options_dict['HELICAL PARAMETERS']['Rise'] = {
                                              'glossory': 'rise',
                                              'default': 0.0,
                                              'validation': lambda x: float(x),
                                              }

_options_dict['HELICAL PARAMETERS']['Buckle'] = {
                                                'glossory': 'buckle',
                                                'default': 0.0,
                                                'validation': lambda x: float(x),
                                                }

_options_dict['HELICAL PARAMETERS']['Propeller'] = {
                                                   'glossory': 'propeller',
                                                   'default': 0.0,
                                                   'validation': lambda x: float(x),
                                                   }

_options_dict['HELICAL PARAMETERS']['Opening'] = {
                                                 'glossory': 'opening',
                                                 'default': 0.0,
                                                 'validation': lambda x: float(x),
                                                 }

_options_dict['HELICAL PARAMETERS']['Shear'] = {
                                               'glossory': 'shear',
                                               'default': 0.0,
                                               'validation': lambda x: float(x),
                                               }

_options_dict['HELICAL PARAMETERS']['Stretch'] = {
                                                 'glossory': 'stretch',
                                                 'default': 0.0,
                                                 'validation': lambda x: float(x),
                                                 }

_options_dict['HELICAL PARAMETERS']['Stagger'] = {
                                                 'glossory': 'stagger',
                                                 'default': 0.0,
                                                 'validation': lambda x: float(x),
                                                 }

_options_dict['HELICAL PARAMETERS']['Inclination'] = {
                                                     'glossory': 'inclination',
                                                     'default': 0.0,
                                                     'validation': lambda x: float(x),
                                                     }

_options_dict['HELICAL PARAMETERS']['Tip'] = {
                                             'glossory': 'tip',
                                             'default': 0.0,
                                             'validation': lambda x: float(x),
                                             }

_options_dict['HELICAL PARAMETERS']['X_Displacement'] = {
                                                        'glossory': 'x-displacement',
                                                        'default': 0.0,
                                                        'validation': lambda x: float(x),
                                                        }

_options_dict['HELICAL PARAMETERS']['Y_Displacement'] = {
                                                        'glossory': 'y-displacement',
                                                        'default': 0.0,
                                                        'validation': lambda x: float(x),
                                                        }

# Runtime Parameters
_options_dict['RUNTIME PARAMETERS'] = {}
_options_dict['RUNTIME PARAMETERS']['Search_Size'] = {
                                                     'glossory': 'number of points to search',
                                                     'default': 36000000,
                                                     'validation': lambda x: int(x),
                                                     }

_options_dict['RUNTIME PARAMETERS']['Force_Field_Type'] = {
                                                          'glossory': 'force field type',
                                                          'default': 'GAFF',
                                                          'validation': lambda x: str(x).upper(),
                                                          }

# Not used yet.
#_options_dict['RUNTIME PARAMETERS']['Force_Field_Parameter_File'] = {
#                                                                    'glossory': 'additional parameter file',
#                                                                    'default': None,
#                                                                    'validation': lambda x: _validate_input_file(x),
#                                                                    }

def _validate_algorithm(algorithm):
    """
    method to validate algorithm. Checks if the requested algorith is available.
    """

    available_algorithms = ['WMC']
    algorithm = str(algorithm).upper()
    if algorithm in available_algorithms:
        return algorithm
    else:
        raise Exception('Available alogrithms are %s' %str(available_algorithms))

_options_dict['RUNTIME PARAMETERS']['Algorithm'] = {
                                                   'glossory': 'three letter code for search algorithm',
                                                   'default': 'WMC',
                                                   'validation': lambda x: _validate_algorithm(x),
                                                   }

_options_dict['RUNTIME PARAMETERS']['Max_Total_Energy'] = {
                                                          'glossory': ('cuttoff for the total energy of the entire molecule' +
                                                                       ' resulting from sequence specified'),
                                                          'default': '1e15',
                                                          'validation': lambda x: float(x),
                                                          }

_options_dict['RUNTIME PARAMETERS']['Max_Angle_Energy'] = {
                                                          'glossory': 'cuttoff for energy from angles between three consecutive atoms',
                                                          'default': '1e15',
                                                          'validation': lambda x: float(x),
                                                          }

_options_dict['RUNTIME PARAMETERS']['Max_Bond_Energy'] = {
                                                         'glossory': 'cuttoff for energy from bond between newly bonded atoms',
                                                         'default': '1e15',
                                                         'validation': lambda x: float(x),
                                                         }

_options_dict['RUNTIME PARAMETERS']['Max_Vdw_Energy'] = {
                                                        'glossory': 'cuttoff for total vdw energy',
                                                        'default': '1e15',
                                                        'validation': lambda x: float(x),
                                                        }

_options_dict['RUNTIME PARAMETERS']['Max_Torsion_Energy'] = {
                                                            'glossory': ('cuttoff for energy from the new torsion formed' +
                                                                         ' bewteen linking backbones'),
                                                            'default': '1e15',
                                                            'validation': lambda x: float(x),
                                                            }

_options_dict['RUNTIME PARAMETERS']['Max_Backbone_Interlink_Distance'] = {
                                                                         'glossory': ('the maximum distance between atom' +
                                                                                      ' linkers in backbone. Conformers that' + 
                                                                                      ' do not meet this criteria are' + 
                                                                                      ' automatically rejected and do not even' + 
                                                                                      ' pass to energy analysis where a chain is built'),
                                                                         'default': 0.15,
                                                                         'validation': lambda x: float(x),
                                                                         }

_options_dict['RUNTIME PARAMETERS']['Strand_Base_Names'] = {
                                                           'glossory': 'list of base names in strand',
                                                           'default': None,
                                                           'validation': lambda x: (tuple([i.strip() for i in x.split(',') if i]) if isinstance(x, str)
                                                                                    else tuple([str(i).strip() for i in tuple(x) if i])),
                                                           }

_options_dict['RUNTIME PARAMETERS']['Is_Double_Stranded'] = {
                                                            'glossory': 'whether to generate double stranded conformers',
                                                            'default': False,
                                                            'validation': lambda x: bool(eval(x.title())) if isinstance(x, str) else bool(x),
                                                            }


for k1 in _options_dict:
    for k2 in _options_dict[k1]:
        _options_dict[k1][k2]['new_value'] = None
        _options_dict[k1][k2]['validated'] = False

# Finished defining default options_dict

def options_dict():
    """
    Function to return a deep copy of options_dict
    """

    return copy.deepcopy(_options_dict) 
