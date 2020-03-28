"""!@file
A file for preparing, explaining and validating options

@namespace options
@brief A file for preparing, explaining and validating options

This option file defines a dictionary containing all the available options for the builder code.
@a options._options_dict has three keys: Backbone, HelicalParameters, and RuntimeParameters.
Each of these contains a dictionary of available options in each category. Each dictionary contains
a short glossory that describe the option, a long glossory giving hints on how to
use this option, a default value, and a validation scheme.

This file also contains additional functions to help in the validation.
"""

from __future__ import division, absolute_import, print_function

import os
from openbabel import openbabel as ob
import yaml

from pnab import __path__


def validate_all_options(options):
    """!@brief A method to validate all options.

    It loops over all available builder options and validates the user input. For options
    that have not been specified by the user, a default value is provided.

    @param options (dict) A dictionary containing all the user-defined options

    @return None; @a options is modified in place

    @sa pNAB.pNAB.__init__
    @sa _options_dict
    """

    # Loop over all available options
    for k1 in _options_dict:
        if k1 not in options:
            if k1 == 'Base': continue
            options[k1] = {}
        for k2 in _options_dict[k1]:
            if k2 not in options[k1]:
                # Add a default value for options not specified by the user
                options[k1][k2] = _options_dict[k1][k2]['default']
            options[k1][k2] = _options_dict[k1][k2]['validation'](options[k1][k2])

    # Validate base options if additional bases are specified
    for k1 in options:
        if 'Base' in k1:
            for k2 in _options_dict['Base']:
                if k2 not in options[k1]:
                     # Add a default value for options not specified by the user
                      options[k1][k2] = _options_dict['Base'][k2]['default']
                options[k1][k2] = _options_dict['Base'][k2]['validation'](options[k1][k2])
            if options[k1].pop('align'):
                _align_nucleobase(options[k1])

    # Validate whether some options are mutually exclusive
    # Cannot build hexad strands for canonical nucleobases
    if options['RuntimeParameters']['is_hexad']:
        for i in ['A', 'G', 'C', 'U', 'T']:
            if i in options['RuntimeParameters']['strand']:
                raise Exception("Cannot build hexads for canonical nucleobases")

def _align_nucleobase(base_options):
    """!@brief Aligns the provided nucleobase to purine or pyrimidine in the nucleic acid base pair standard reference frame

    This function aligns the provided nucleobase to either guanine if it has two rings or to uracil if it has one ring.
    The provided linker atoms that connect the base to the backbone are used for alignment. A third atom used in the alignment is 
    determined by rearraning the atoms in the molecules using the canonical order. The function loops over all atoms in the provided
    nucleobase and the reference nucleobase and if two atoms have the same atomic number then these atoms are used for the alignment.
    A new file with the aligned nucleobase is created and used instead of the provided file.

    @param base_options (dict) The nucleobase options defined in the input file
    """

    conv = ob.OBConversion()
    nucleobase = ob.OBMol()
    # Read the provided nucleobase
    conv.ReadFile(nucleobase, base_options['file_path'])
    rings = list(nucleobase.GetSSSR())
    # Use either guanine or uracil as a reference
    if len(rings) == 2:
        reference = 'G'
    elif len(rings) == 1:
        reference = 'U'
    else:
        return

    # Read the data for the reference nucleobase
    ref_library = yaml.load(open(os.path.join(__path__[0], 'data', 'bases_library.yaml')), yaml.FullLoader) 
    ref_linker = ref_library["Base " + reference]['linker']
    ref_path = ref_library["Base " + reference]['file_path']

    ref = ob.OBMol()
    conv.ReadFile(ref, os.path.join(__path__[0], 'data', ref_path))

    # Find the ring atoms in the reference base
    # compute the centroid
    centroid_ref = ob.vector3()
    num_atoms = 0
    ref_ring = ob.OBMol()
    # Add the two linker atoms
    ref_ring.AddAtom(ref.GetAtom(ref_linker[1]))
    ref_ring.AddAtom(ref.GetAtom(ref_linker[0]))
    ref_ring.AddBond(1, 2, 1)
    for atom in ob.OBMolAtomIter(ref):
        if atom.IsInRing() and atom.GetIdx() != ref_linker[0]:
            ref_ring.AddAtom(atom)
            centroid_ref += atom.GetVector()
            num_atoms += 1
    centroid_ref /= num_atoms

    # Find the ring atoms in the provided nucleobase
    # compute the centroid
    centroid_nucleobase = ob.vector3()
    num_atoms = 0 
    nucleobase_ring = ob.OBMol()
    # Add the two linker atoms
    nucleobase_ring.AddAtom(nucleobase.GetAtom(base_options['linker'][1]))
    nucleobase_ring.AddAtom(nucleobase.GetAtom(base_options['linker'][0]))
    nucleobase_ring.AddBond(1, 2, 1)
    for atom in ob.OBMolAtomIter(nucleobase):
        if atom.IsInRing() and atom.GetIdx() != base_options['linker'][0]:
            nucleobase_ring.AddAtom(atom)
            centroid_nucleobase += atom.GetVector()
            num_atoms += 1
    centroid_nucleobase /= -1*num_atoms

    # Pick three atoms for alignment
    # First, add the two linker atoms
    ref_align = ob.OBMol()
    nucleobase_align = ob.OBMol()
    for i in range(1, 3):
        ref_align.AddAtom(ref_ring.GetAtom(i))
        nucleobase_align.AddAtom(nucleobase_ring.GetAtom(i))

    # Reorder the atoms using the canonical order
    canonical = ob.OBOp.FindType("canonical")
    canonical.Do(ref_ring)
    canonical.Do(nucleobase_ring)

    # Delete the two linker atoms
    for atom in ob.OBAtomAtomIter(ref_ring.GetAtom(1)):
        ref_ring.DeleteAtom(atom)
    ref_ring.DeleteAtom(ref_ring.GetAtom(1))

    for atom in ob.OBAtomAtomIter(nucleobase_ring.GetAtom(1)):
        nucleobase_ring.DeleteAtom(atom)
    nucleobase_ring.DeleteAtom(nucleobase_ring.GetAtom(1))

    # Search for two atoms that have the same index and atomic number
    for i in range(1, ref_ring.NumAtoms() + 1):
        if ref_ring.GetAtom(i).GetAtomicNum() == nucleobase_ring.GetAtom(i).GetAtomicNum():
            ref_align.AddAtom(ref_ring.GetAtom(i))
            nucleobase_align.AddAtom(nucleobase_ring.GetAtom(i))
            break

    # Align the three atoms in the reference nucleobase and the provided nucleobase
    align = ob.OBAlign(ref_align, nucleobase_align, True, True)
    align.Align()

    # Apply the rotation to all the atoms
    matrix = align.GetRotMatrix()
    array = ob.doubleArray(9)
    matrix.GetArray(array)
    nucleobase.Translate(centroid_nucleobase)
    nucleobase.Rotate(array)
    nucleobase.Translate(centroid_ref)

    # Write a new file with the aligned coordinates and update the options
    base_options['file_path'] = os.path.splitext(base_options['file_path'])[0] + '_aligned' + os.path.splitext(base_options['file_path'])[1]
    conv.WriteFile(nucleobase, base_options['file_path'])

def _validate_input_file(file_name):
    """!@brief Method to validate that the given file exists.

    Check whether the file exists or not. If it does not exist as it is, the function checks
    whether the file is in the "pnab/data" directory. 

    @param file_name (str) Path to a file

    @return file_name after validation

    @sa validate_all_options
    @sa _options_dict
    """

    file_name = str(file_name)
    if not os.path.isfile(file_name):
        # Check if the file exists in the "pnab/data" directory
        file_name2 = os.path.join(__path__[0], 'data', file_name)
        if not os.path.isfile(file_name2):
            raise Exception("Cannot find file: %s" %file_name)
        file_name = file_name2

    return file_name


def _validate_atom_indices(indices):
    """!@brief Method to validate provided lists of atom indices.

    Check whether the list has two atom indices

    @param indices (str|list) atom indices

    @return indices after validation

    @sa validate_all_options
    @sa _options_dict
    """

    indices = list(eval(indices)) if isinstance(indices, str) else list(indices)
    if len(indices) != 2:
        raise Exception("Incorrect number of atoms. Must provide two indices.")

    for i in range(len(indices)):
        indices[i] = int(indices[i])
        if indices[i] == 0:
            raise Exception("Use 1-based index")

    return indices


def _validate_helical_parameters(hp_i):
    """!@brief Method to validate provided helical parameters for each parameter.

    Check whether the list has correct helical parameter specifications.
    The correct specifications are [initial value in a range, final value in a range, number of steps].
    If a single value is provided, then we assume it is the only value in a range.
    if two values are provided, then we assume it is a range with one configuration.

    This specification is for internal use in the driver for generating multiple configurations
    and running them in parallel. The C++ code accepts a single value for each helical parameter.

    @param  hp_i (str|float|list) Specifications for helical parameter i

    @return hp_i after validation

    @sa validate_all_options
    @sa _options_dict
    @sa pNAB.pNAB.run
    """

    if isinstance(hp_i, str):
        hp_i = eval(hp_i)

    # If it is a single int or float value
    if isinstance(hp_i, (float, int)):
        hp_i = [hp_i, hp_i, 1]
    # if it is a list
    else:
        hp_i = list(hp_i)

    # If a list has one value, then we assume it is a single value for the helical parameters
    if len(hp_i) == 1:
        hp_i = [hp_i[0], hp_i[0], 1]
    # If the list has two values, then we assume it is a range with one generated helical configuration
    elif len(hp_i) == 2:
        # Default of 1 value tested for a range
        hp_i.append(1)
    elif len(hp_i) > 3:
        raise Exception("Helical parameters must have at most three values")

    return hp_i


def _validate_energy_filter(energy_filter):
    """!@brief Method to validate provided energy filter.

    Energy filter has five thresholds [bond, angle, torsion, van der Waals, total]

    @param energy_filter (str|list) A list of five energy thresholds

    @return energy_filter after validation

    @sa validate_all_options
    @sa _options_dict
    """

    if isinstance(energy_filter, str):
        energy_filter = eval(energy_filter)

    energy_filter = list(energy_filter)
    if len(energy_filter) != 5:
        raise Exception("Five values must be provided for the energy filter")

    for i in energy_filter:
        if not isinstance(i, (float, int)):
            raise Exception("Provide valid numbers")

    return energy_filter

def _validate_strand(strand):
    """!@brief Method to validate provided strand sequence.

    We use FASTA strings to specify the sequence in the strand.
    This is for use in the driver. The C++ code accepts base
    names with more than one letter. However, writing a FASTA
    sequence is easier.

    The names of the bases are converted to upper case letters

    @param strand (str|list) FASTA sequence of the strand

    @return strand list after validation

    @sa validate_all_options
    @sa _options_dict
    """

    strand = list(strand)
    for i in range(len(strand)):
        strand[i] = strand[i].upper()
    # Check if both canonical and noncanonical bases are in the sequence
    if "M" in strand or "Y" in strand:
        for i in ["A", "G", "C", "T", "U"]:
            if i in strand:
                raise Exception("Cannot combine canonincal and non-canonical nucleobases")

    return strand

def _validate_bool_list(x):
    """!@brief Method to validate a list of bools

    @param x (list) list of bools

    @return validated list of bools

    @sa validate_all_options
    @sa _options_dict
    """

    x = list(x)
    for i, val in enumerate(x):
        x[i] = bool(eval(val.title())) if isinstance(val, str) else bool(val)

    return x


# Set glossory of options, default values and validation methods

## @brief Nucleic acid builder options used in the driver.
#
# This is an internal options dictionary to be used for listing, explaining
# and validating the options used in the code. All options here have corresponding
# options in the C++ code, except @a options._options_dict['RuntimeParameters']['pair_a_u'],
# which is used in the driver to override the default adenine-thymine pairing.
# @sa jupyter_widgets.builder()
# @sa PNAB::Backbone
# @sa PNAB::HelicalParameters
# @sa PNAB::RuntimeParameters
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
                                         'default': '',
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
                                           'validation': lambda x: [_validate_atom_indices(i) for i in list(x)],
                                           }

# Base Parameters
_options_dict['Base'] = {}

_options_dict['Base']['file_path'] = {
                                       'glossory': 'Path to file containing the molecular structure of the base',
                                       'long_glossory': ('Path to the file containing the three-dimensional' +
                                                        ' molecular structure of the nucleobase (e.g. PDB file).' +
                                                        ' The base must contain hydrogen atoms and not contain' +
                                                        ' atoms from the backbone. The base must be in the correct standard frame of reference.'),
                                       'default': '',
                                       'validation': lambda x: _validate_input_file(x),
                                       }
_options_dict['Base']['linker'] = {
                                    'glossory': 'Two atoms forming a vector connecting to backbone',
                                    'long_glossory': ('Select two atoms that form the vector connecting the nucleobase molecule to' +
                                                      ' the backbone. The terminal atom will be deleted as the bond with the nucleobase' +
                                                      ' is formed.'),
                                    'default': [0, 0],
                                    'validation': lambda x: _validate_atom_indices(x),
                                    }
_options_dict['Base']['code'] = {
                                  'glossory': 'Three-letter code',
                                  'long_glossory': 'This code is used as the residue name in the output PDB files.',
                                  'default': 'RES',
                                  'validation': lambda x: str(x),
                                  }
_options_dict['Base']['name'] = {
                                  'glossory': 'One-letter base name',
                                  'long_glossory': ('This name is used when specifying the strand sequence. It must not be one of the' +
                                                    ' names defined in the program library (A, G, C, T, U, M, Y).'),
                                  'default': 'R',
                                  'validation': lambda x: str(x),
                                  }
_options_dict['Base']['pair_name'] = {
                                       'glossory': 'Name of the pairing base',
                                       'long_glossory': 'The one-letter name of the pairing base.',
                                       'default': '',
                                       'validation': lambda x: str(x),
                                       }
_options_dict['Base']['align'] = {
                                 'glossory': 'Align the base to the standard frame of reference',
                                 'long_glossory': ('Aligns the provided nucleobase to either purine or pyrimidine in the base pair standard frame of reference.'),
                                 'default': False,
                                 'validation': lambda x: bool(x),
                                 }

# Helical Parameters
_options_dict['HelicalParameters'] = {}
_options_dict['HelicalParameters']['h_twist'] = {
                                               'glossory': 'Helical Twist (degree)',
                                               'long_glossory': ('Select a single value or a range of values. If you select a range of values,' + 
                                                                 ' you can control the number of configurations generated in that range by increasing' + 
                                                                 ' the number of steps. Uniformly spaced configurations in the provided range will be generated.'), 
                                               'default': [0.0, 0.0, 1],
                                               'validation': lambda x: _validate_helical_parameters(x),
                                               }
_options_dict['HelicalParameters']['inclination'] = {
                                                    'glossory': 'Inclination (degree)',
                                                    'long_glossory': ('Select a single value or a range of values. If you select a range of values,' + 
                                                                      ' you can control the number of configurations generated in that range by increasing' + 
                                                                      ' the number of steps. Uniformly spaced configurations in the provided range will be generated.' +
                                                                      ' Not defined for the hexad geometry.'), 
                                                    'default': [0.0, 0.0, 1],
                                                    'validation': lambda x: _validate_helical_parameters(x),
                                                    }
_options_dict['HelicalParameters']['tip'] = {
                                            'glossory': 'Tip (degree)',
                                            'long_glossory': ('Select a single value or a range of values. If you select a range of values,' + 
                                                              ' you can control the number of configurations generated in that range by increasing' + 
                                                              ' the number of steps. Uniformly spaced configurations in the provided range will be generated.' + 
                                                              ' Not defined for the hexad geometry.'), 
                                            'default': [0.0, 0.0, 1],
                                            'validation': lambda x: _validate_helical_parameters(x),
                                            }
_options_dict['HelicalParameters']['h_rise'] = {
                                                'glossory': 'Helical Rise (Angstrom)',
                                                'long_glossory': ('Select a single value or a range of values. If you select a range of values,' + 
                                                                  ' you can control the number of configurations generated in that range by increasing' + 
                                                                  ' the number of steps. Uniformly spaced configurations in the provided range will be generated.'), 
                                                'default': [0.0, 0.0, 1],
                                                'validation': lambda x: _validate_helical_parameters(x),
                                                }
_options_dict['HelicalParameters']['x_displacement'] = {
                                                       'glossory': 'X-Displacement (Angstrom)',
                                                       'long_glossory': ('Select a single value or a range of values. If you select a range of values,' + 
                                                                         ' you can control the number of configurations generated in that range by increasing' + 
                                                                         ' the number of steps. Uniformly spaced in the provided range will be generated.' + 
                                                                         ' Not defined for the hexad geometry.'), 
                                                       'default': [0.0, 0.0, 1],
                                                       'validation': lambda x: _validate_helical_parameters(x),
                                                       }
_options_dict['HelicalParameters']['y_displacement'] = {
                                                       'glossory': 'Y-Displacement (Angstrom)',
                                                       'long_glossory': ('Select a single value or a range of values. If you select a range of values,' + 
                                                                         ' you can control the number of configurations generated in that range by increasing' + 
                                                                         ' the number of steps. Uniformly spaced configurations in the provided range will be generated.' + 
                                                                         ' Not defined for the hexad geometry.'), 
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
                                                  'glossory': 'Number of iterations to search over dihedral angles',
                                                  'long_glossory': ('Number of iterations (or generations in genetic algorithm) for' +
                                                                    ' searching over dihedral angles. This should be chosen such that it samples' +
                                                                    ' the dihedral angles enough while keeping the time of the search reasonbly short.' +
                                                                    ' This should be increased for systems with greater number of rotatable bonds in the backbone.'),
                                                  'default': 1000000,
                                                  'validation': lambda x: int(x),
                                                  }
_options_dict['RuntimeParameters']['seed'] = {
                                             'glossory': 'The seed for the random number generator',
                                             'long_glossory': 'Use the same value for the seed to get reproducible results in the same computer platform.',
                                             'default': 0,
                                             'validation': lambda x: int(abs(x)),
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
                                                              'default': 300.0,
                                                              'validation': lambda x: float(x),
                                                              }
_options_dict['RuntimeParameters']['monte_carlo_temperature'] = {
                                                                'glossory': 'Temperature used in the Monte Carlo procedure (K)',
                                                                'long_glossory': ('This temperature controls the acceptance and rejection ratio of the Monte Carlo steps'),
                                                                'default': 300.0,
                                                                'validation': lambda x: float(x),
                                                                }
_options_dict['RuntimeParameters']['population_size'] = {
                                                        'glossory': 'The size of the population',
                                                        'long_glossory': ('The size of the population in the genetic algorithm search.'),
                                                        'default': 1000,
                                                        'validation': lambda x: int(x),
                                                        }
_options_dict['RuntimeParameters']['mutation_rate'] = {
                                                      'glossory': 'Mutation rate',
                                                      'long_glossory': ('Mutation rate in the genetic algorithm search. Used to intorduce new values' +
                                                                        ' for the dihedral angles in the population.'),
                                                      'default': 0.5,
                                                      'validation': lambda x: float(x),
                                                      }
_options_dict['RuntimeParameters']['crossover_rate'] = {
                                                       'glossory': 'Crossover rate',
                                                       'long_glossory': ('Crossover or mating rate in the genetic algorithm search.' +
                                                                         ' Used to exchange dihedral angles between individuals in the population.'),
                                                       'default': 0.5,
                                                       'validation': lambda x: float(x),
                                                       }
_options_dict['RuntimeParameters']['ff_type'] = {
                                                'glossory': 'Force field type',
                                                'long_glossory': 'Force field for computing the energy of the system.', 
                                                'default': 'MMFF94',
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
                                                     'default': 0.2,
                                                     'validation': lambda x: float(x),
                                                     }
_options_dict['RuntimeParameters']['energy_filter'] = {
                                                      'glossory': ('Maximum energy for newly formed bonds in the backbone (kcal/mol)\n' +
                                                                   'Maximum energy for newly formed angles in the backbone (kcal/mol)\n' +
                                                                   'Maximum torsional energy for rotatable bonds (kcal/mol/nucleotide)\n' +
                                                                   'Maximum van der Waals energy (kcal/mol/nucleotide)\n' +
                                                                   'Maximum total energy (kcal/mol/nucleotide)\n'),
                                                      'long_glossory': ('This is the bond stretching energy for the newly formed bonds between the first two adjacent backbone molecules.\n' +
                                                                        'This is the angle bending energy for the newly formed angles between the first two adjacent backbone molecules.\n' +
                                                                        'This is the torsional energy for all the rotatable backbone bonds.\n' + 
                                                                        'This is the total van der Waals energy of the system.\n' + 
                                                                        'This is the total energy of the system.\n'),
                                                      'default': (1, 4, 10, 500, 10000000000),
                                                      'validation': lambda x: _validate_energy_filter(x), 
                                                      }
_options_dict['RuntimeParameters']['strand'] = {
                                               'glossory': 'FASTA string for nucleotide sequence (e.g. GCAT or MYMY) ',
                                               'long_glossory': ('See the Bases section for the defined bases and their one-letter abbreviation.' + 
                                                                 'The canonical nucleobases cannot be mixed with the nucleobases that form the hexad geometries' +
                                                                 ' because they have different standard frame of reference. The names are case-insensitive.'),
                                               'default': None,
                                               'validation': lambda x: _validate_strand(x),
                                               }
_options_dict['RuntimeParameters']['pair_A_U'] = {
                                                 'glossory': 'Pair A with U (Default is A-T pairing)',
                                                 'long_glossory': 'Check box if you want to pair adenine with uracil instead of thymine.',
                                                 'default': False,
                                                 'validation': lambda x: bool(eval(x.title())) if isinstance(x, str) else bool(x),
                                                 }
_options_dict['RuntimeParameters']['is_hexad'] = {
                                                 'glossory': 'Hexad strands',
                                                 'long_glossory': 'Check box if you want to apply the 60 degrees rotation for the hexad strands.',
                                                 'default': False,
                                                 'validation': lambda x: bool(eval(x.title())) if isinstance(x, str) else bool(x),
                                                 }
_options_dict['RuntimeParameters']['build_strand'] = {
                                                     'glossory': 'Build the strand',
                                                     'long_glossory': ('Check the first box only if you want to build a single strand.' + 
                                                                       ' Check the first two boxes if you want to build a duplex.' + 
                                                                       ' Check all boxes if you want to build a hexad. Canonical nucleobases' + 
                                                                       ' can only build a single strand or a duplex.'),
                                                     'default': [True, False, False, False, False, False],
                                                     'validation': lambda x: _validate_bool_list(x),
                                                     }
_options_dict['RuntimeParameters']['strand_orientation'] = {
                                                           'glossory': 'Orientation of each strand (up or down)',
                                                           'long_glossory': ('Orientation of each strand. For DNA and RNA, check the first box' +
                                                                             ' and uncheck the second box to build anti-parallel duplex. For hexads' + 
                                                                             ' Check the boxes to build any combination of parallel and anti-parallel strands.'),
                                                           'default': [True, True, True, True, True, True],
                                                           'validation': lambda x: _validate_bool_list(x),
                                                           }
_options_dict['RuntimeParameters']['glycosidic_bond_distance'] = {
                                                                 'glossory': 'The distance of the glycosidic bond or its equivalent',
                                                                 'long_glossory': ("A user-defined distance for the glycosidic bond" + 
                                                                                   " If it is zero, the van der Waals radii are used to" + 
                                                                                   " determine an appropriate distance. This distance is set for" + 
                                                                                   " first nucleotide. The distances in the other nucleotides are" + 
                                                                                   " set by the periodic condition."),
                                                                 'default': 0.0,
                                                                 'validation': lambda x: float(x),
                                                                 }
_options_dict['RuntimeParameters']['num_candidates'] = {
                                                       'glossory': 'Quit after finding the specified number of accepted candidates',
                                                       'long_glossory': ('Quit after finding the specified number of accepted candidates for a given helical configuration.' + 
                                                                         ' Can be used to save time if only a few backbone conformations are needed.' +
                                                                         ' If the specified number of candidates is not satisfied, the search will continue' + 
                                                                         ' until the requested number of steps is completed.'),
                                                       'default': 10,
                                                       'validation': lambda x: int(x),
                                                       }
