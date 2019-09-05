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
    if len(x) != 5:
        raise Exception("Five values must be provided for the energy filter")

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

def _validate_fixed_bonds(x):
    x = list(x)
    for i in range(len(x)):
        x[i] = list(x[i])
        if len(x[i]) != 2:
            raise Exception("Specify two atoms for fixed rotatable bonds")

    return x


# Set glossory of options, default values and validation methods
_options_dict = {}

# Backbone Parameter
_options_dict['Backbone'] = {}
_options_dict['Backbone']['file_path'] = {
                                         'glossory': ('Path to the file containing the molecular' +
                                                     ' structure of the backbone'),
                                         'long_glossory': ('Path to the file containing the three-dimensional' + 
                                                          ' molecular structure of the backbone (e.g. PDB file).' + 
                                                          ' The backbone must contain hydrogen atoms and not contain' +
                                                          ' atoms from the nucleobase. The atoms included here must' + 
                                                          ' be for one nucleotide.'),
                                         'default': 'backbone.pdb',
                                         'validation': lambda x: _validate_input_file(x),
                                         }
_options_dict['Backbone']['interconnects'] = {
                                             'glossory': 'Two atoms connecting to the two backbones',
                                             'long_glossory': ('Select two atoms that connect the backbone molecule to' +
                                                               ' the two neighboring backbone molecules. Extra hydrogen atoms' + 
                                                               ' connected to these two linker atoms are deleted. The order of' + 
                                                               ' atoms control the directionality of the backbone.'),
                                             'default': None,
                                             'validation': lambda x: _validate_atom_indices(x),
                                             }
_options_dict['Backbone']['linker'] = {
                                      'glossory': ('Two atoms forming the vector' + 
                                                  ' connecting to base'),
                                      'long_glossory': ('Select two atoms that form the vector connecting the backbone molecule to' +
                                                        ' the nucleobase. The terminal atom will be deleted as the bond with the nucleobase' +
                                                        ' is formed.'),
                                      'default': None,
                                      'validation': lambda x: _validate_atom_indices(x),
                                      }
_options_dict['Backbone']['fixed_bonds'] = {
                                           'glossory': ('Indices of fixed rotatable dihedral bonds'),
                                           'long_glossory': ('Select the indices of the two atoms that are at the center' + 
                                                             ' of the rotatable dihedral angle.'),
                                           'default': [],
                                           'validation': lambda x: _validate_fixed_bonds(x),
                                           }

# Helical Parameters
_options_dict['HelicalParameters'] = {}
_options_dict['HelicalParameters']['h_twist'] = {
                                               'glossory': 'Helical Twist (degree)',
                                               'long_glossory': ('Select a single value or a range of values. If you select a range of values,' + 
                                                                 ' you can control the number of configurations generated in that range by increasing' + 
                                                                 ' the number of steps. Random configurations in the provided range will be generated.'), 
                                               'default': [0.0, 0.0, 1],
                                               'validation': lambda x: _validate_helical_parameters(x),
                                               }
_options_dict['HelicalParameters']['inclination'] = {
                                                    'glossory': 'Inclination (degree)',
                                                    'long_glossory': ('Select a single value or a range of values. If you select a range of values,' + 
                                                                      ' you can control the number of configurations generated in that range by increasing' + 
                                                                      ' the number of steps. Random configurations in the provided range will be generated.'), 
                                                    'default': [0.0, 0.0, 1],
                                                    'validation': lambda x: _validate_helical_parameters(x),
                                                    }
_options_dict['HelicalParameters']['tip'] = {
                                            'glossory': 'Tip (degree)',
                                            'long_glossory': ('Select a single value or a range of values. If you select a range of values,' + 
                                                              ' you can control the number of configurations generated in that range by increasing' + 
                                                              ' the number of steps. Random configurations in the provided range will be generated.'), 
                                            'default': [0.0, 0.0, 1],
                                            'validation': lambda x: _validate_helical_parameters(x),
                                            }
_options_dict['HelicalParameters']['h_rise'] = {
                                                'glossory': 'Helical Rise (Angstrom)',
                                                'long_glossory': ('Select a single value or a range of values. If you select a range of values,' + 
                                                                  ' you can control the number of configurations generated in that range by increasing' + 
                                                                  ' the number of steps. Random configurations in the provided range will be generated.'), 
                                                'default': [0.0, 0.0, 1],
                                                'validation': lambda x: _validate_helical_parameters(x),
                                                }
_options_dict['HelicalParameters']['x_displacement'] = {
                                                       'glossory': 'X-Displacement (Angstrom)',
                                                       'long_glossory': ('Select a single value or a range of values. If you select a range of values,' + 
                                                                         ' you can control the number of configurations generated in that range by increasing' + 
                                                                         ' the number of steps. Random configurations in the provided range will be generated.'), 
                                                       'default': [0.0, 0.0, 1],
                                                       'validation': lambda x: _validate_helical_parameters(x),
                                                       }
_options_dict['HelicalParameters']['y_displacement'] = {
                                                       'glossory': 'Y-Displacement (Angstrom)',
                                                       'long_glossory': ('Select a single value or a range of values. If you select a range of values,' + 
                                                                         ' you can control the number of configurations generated in that range by increasing' + 
                                                                         ' the number of steps. Random configurations in the provided range will be generated.'), 
                                                       'default': [0.0, 0.0, 1],
                                                       'validation': lambda x: _validate_helical_parameters(x),
                                                       }

# Runtime Parameters
_options_dict['RuntimeParameters'] = {}
_options_dict['RuntimeParameters']['search_algorithm'] = {
                                                         'glossory': 'Search algorithm',
                                                         'long_glossory': ('There are six search algorithms comprising four classes:\n' + 
                                                                           '1) Monte Carlo search, 2) Random Search, 3) Genetic Algorithm Search, and 4) Systematic Search\n'
                                                                           'The first three algorithms are not deterministic. Weighted algorithms do not use' + 
                                                                           ' uniform distributions for generating random dihedral angles. Instead, the probability' +
                                                                           ' distribution for each dihedral angle is weighted by exp(-Ei/kT)/sum(exp(-E/kT) where' + 
                                                                           ' Ei is the torsional energy at dihedral angle i.'
                                                                          ),
                                                         'default': 'weighted monte carlo search',
                                                         'validation': lambda x: x.lower(),
                                                         }
_options_dict['RuntimeParameters']['num_steps'] = {
                                                  'glossory': 'Number of points to search over dihedral angles',
                                                  'long_glossory': ('Number of iterations (or generations in genetic algorithm) for' +
                                                                    ' searching over dihedral angles. This should be chosen such that it samples' +
                                                                    ' the dihedral angles enough while keeping the time of the search reasonbly short.' +
                                                                    ' This should be increased for systems with greater number of rotatable bonds in the backbone.'),
                                                  'default': 1000000,
                                                  'validation': lambda x: int(x),
                                                  }
_options_dict['RuntimeParameters']['dihedral_step'] = {
                                                      'glossory': 'Dihedral step size for systematic search (degree)',
                                                      'long_glossory': ('The number of search points for systematic search is given by:'
                                                                        ' (360/step)^(number of rotatable bonds).'),
                                                      'default': 2,
                                                      'validation': lambda x: float(x),
                                                      }
_options_dict['RuntimeParameters']['weighting_temperature'] = {
                                                              'glossory': 'Weighting temperature (K)',
                                                              'long_glossory': ('Temperature used for weighting the probability of each dihedral angle'),
                                                              'default': 298.0,
                                                              'validation': lambda x: float(x),
                                                              }
_options_dict['RuntimeParameters']['monte_carlo_temperature'] = {
                                                                'glossory': 'Temperature used in the Monte Carlo procedure (K)',
                                                                'long_glossory': ('This temperature controls the acceptance and rejection ratio of the Monte Carlo steps'),
                                                                'default': 298.0,
                                                                'validation': lambda x: float(x),
                                                                }
_options_dict['RuntimeParameters']['type'] = {
                                             'glossory': 'Force field type',
                                             'long_glossory': 'Force field for computing the energy of the system.', 
                                             'default': 'GAFF',
                                             'validation': lambda x: str(x).upper(),
                                             }
_options_dict['RuntimeParameters']['max_distance'] = {
                                                     'glossory': ('The maximum distance between atom' +
                                                                  ' linkers in backbone (Angstrom)' 
                                                                  ),
                                                     'long_glossory': ('This is the distance between the second terminal atom in one backbone' + 
                                                                       ' and the first terminal atom in the adjacent backbone. When the distance' + 
                                                                       ' between these two atoms is below this threshold, a bond can form between' + 
                                                                       ' the two backbone molecules. This should be below 0.1 Angstroms. Conformers' +
                                                                       ' that do not pass this threshold are immediately rejected without proceeding to build' +
                                                                       ' the system.'),
                                                     'default': 0.1,
                                                     'validation': lambda x: float(x),
                                                     }
_options_dict['RuntimeParameters']['energy_filter'] = {
                                                      'glossory': ('Maximum energy for newly formed bonds in the backbone (kcal/mol/bond)\n' +
                                                                   'Maximum energy for newly formed angles in the backbone (kcal/mol/angle)\n' +
                                                                   'Maximum torsional energy for rotatable bonds (kcal/mol/nucleotide)\n' +
                                                                   'Maximum van der Waals energy (kcal/mol/nucleotide)\n' +
                                                                   'Maximum total energy (kcal/mol/nucleotide)\n'),
                                                      'long_glossory': ('This is the bond stretching energy for the newly formed bonds between two adjacent backbone molecules.\n' +
                                                                        'This is the angle bending energy for the newly formed angles between two adjacent backbone molecules.\n' +
                                                                        'This is the torsional energy for all the rotatable backbone bonds.\n' + 
                                                                        'This is the total van der Waals energy of the system.\n' + 
                                                                        'This is the total energy of the system.\n'),
                                                      'default': (2, 2, 5, 0, 10000000000),
                                                      'validation': lambda x: _validate_energy_filter(x), 
                                                      }
_options_dict['RuntimeParameters']['strand'] = {
                                               'glossory': 'FASTA string for nucleotide sequence (e.g. GCAT or XYXY) ',
                                               'long_glossory': ('See the Bases section for the defined bases and their one-letter abbreviation.' + 
                                                                 'The canonical nucleobases cannot be mixed with the nucleobases that form the hexad geometries' +
                                                                 ' because they have different standard frame of reference. The names are case-insensitive.'),
                                               'default': None,
                                               'validation': lambda x: _validate_strand(x),
                                               }
_options_dict['RuntimeParameters']['is_double_stranded'] = {
                                                           'glossory': 'Double strands',
                                                           'long_glossory': 'Check box if you want to build double strands for the canonical nucleobases.', 
                                                           'default': False,
                                                           'validation': lambda x: bool(eval(x.title())) if isinstance(x, str) else bool(x),
                                                           }
_options_dict['RuntimeParameters']['pair_A_U'] = {
                                                 'glossory': 'Pair A with U (Default is A-T pairing)',
                                                 'long_glossory': 'Check box if you want to pair adenine with uracil instead of thymine.',
                                                 'default': False,
                                                 'validation': lambda x: bool(eval(x.title())) if isinstance(x, str) else bool(x),
                                                 }
_options_dict['RuntimeParameters']['is_hexad'] = {
                                                 'glossory': 'Hexad strands',
                                                 'long_glossory': 'Check box if you want to build hexad strands for the noncanonical nucleobases.',
                                                 'default': False,
                                                 'validation': lambda x: bool(eval(x.title())) if isinstance(x, str) else bool(x),
                                                 }
_options_dict['RuntimeParameters']['strand_orientation'] = {
                                                           'glossory': 'Orientation of each strand in the hexad (up or down)',
                                                           'long_glossory': ('Orientation of the strands in the hexad. This can be used to build' +
                                                                            '  parallel or antiparallel configurations.'),
                                                           'default': [True, True, True, True, True, True],
                                                           'validation': lambda x: _validate_strand_orientation(x),
                                                           }
