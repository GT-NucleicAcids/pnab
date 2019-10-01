"""!@file
A file for displaying widgets in the Jupyter notebook

@namespace jupyter_widgets
@brief A file for displaying widgets in the Jupyter notebook

This file containts widgets for sepecifying options in the Jupyter
notebook using the ipywidgets library (https://github.com/jupyter-widgets/ipywidgets).
Using Jupyter notebook is not necessary for using the program, as
the program can be run in as a python script. However, the widgets
provide a graphical user interface for making input files, running the code,
and displaying the results. The Jupyter notebook can be run locally
or on the clouds using Binder (https://mybinder.org/v2/gh/alenaizan/pnab/master?filepath=binder).

The widgets have four main components: Backbone, Bases, Helical Parameters,
and Runtime Parameters. The Bases section does not require the user input
but displays some information about the available nucleobases in the library.
The NGLView library is used for visualization

To run the widgets from a Jupyter notebook, simply execute the following:

@code
import pnab
pnab.builder()
@endcode
"""

import os
import copy
import datetime
from zipfile import ZipFile

import yaml
import ipywidgets as widgets
from IPython.display import display, Javascript, Image

from pnab.driver.driver import pNAB
from pnab.driver.options import _options_dict
from pnab import __path__

##@brief Stores the widgets corresponding to the user-defined options
#
# This dictionary has the same structure as the options._options_dict
# dictionary that has all the available options for the code. 
# The widgets are added to this dictionary and the values of
# the widgets are extracted when the user runs the code
#
# @sa options._options_dict
# @sa extract_options

input_options = {}

def view_nglview(molecule, label=False):
    """!@brief Display molecules using the NGLView viewer.

    A function to view molecules using the NGLView project (https://github.com/arose/nglview).
    It is used to display the geometry of the backbone with atom index labels. It is
    also used to display accepted nucleic acid candidates generated after
    the search. For backbone molecules, the number of atoms in the backbone
    is computed here and returned to be used for bounding the displayed atomic
    indices in backbone widgets.


    @param molecule (str) Path to a file containing the 3D structure of the molecule
    @param label (bool) Whether to display atom index labels

    @return number of atoms in @a molecule if @a label is True, else None

    @attention If the molecules are not displayed correctly, you may need to execute
        @code{.sh} jupyter-nbextension enable nglview --py --sys-prefix @endcode

    @sa jupyter_widgets.path
    @sa jupyter_widgets.single_result
    """

    import openbabel
    import nglview

    # Check if the file exists
    if not os.path.isfile(molecule):
        return 0

    # If a label is requested, then this is a backbone molecule
    if label:
        # Read the backbone molecule and convert it to a pdb format
        mol = openbabel.OBMol()
        conv = openbabel.OBConversion()
        fmt = conv.FormatFromExt(molecule)
        conv.SetInAndOutFormats(fmt.GetID(), 'pdb')
        conv.ReadFile(mol, molecule)

        # Extract the number of atoms
        num_atoms = mol.NumAtoms()
        mol = conv.WriteString(mol)

        # Display the backbone molecule
        struct = nglview.TextStructure(mol)
        view = nglview.NGLWidget(struct, defaultRepresentation=False)
        view.camera = 'orthographic'
        view.add_representation('licorice')
        view.add_representation('label', labelType='serial', backgroundColor='black', showBackground=True)
        view.center()
        display(view)

        # Return the number of atoms to be used for the backbone widgets
        return num_atoms

    # Show accepted candidates
    else:
        struct = nglview.FileStructure(molecule)
        view = nglview.NGLWidget(struct, defaultRepresentation=False)
        view.camera = 'orthographic'
        view.add_representation('licorice')
        view.center()
        display(view)


def fixed_bonds(num_bonds, num_atoms, param):
    """!@brief Display widgets for determining indices of fixed bonds

    This function is changed interactively bases on the number of fixed
    bonds requested by the user.

    @param num_bonds (int) number of fixed bonds; changed dynamically by the user
    @param num_atoms (int) number of atoms in the backbone
    @param param (dict) @a options._options_dict['Backbone']

    @returns None; @a jupyter_widgets.input_options is modified in place

    @sa path 
    @sa input_options
    @sa options._options_dict
    """

    # Add widget to the widget dictionary
    input_options['Backbone']['fixed_bonds'] = []

    # Capture output. This is used to delete widget when no fixed bond is requested
    out = widgets.Output()
    display(out)
    if num_bonds == 0:
        out.clear_output()

    # Loop over the number of requested fixed bonds
    for i in range(num_bonds):
        try:
            # Check to see if values exist in the used input file
            v1, v2 = param['fixed_bonds']['default'][i][0], param['fixed_bonds']['default'][i][1]
        except IndexError:
            # if they don't exist, then set them to 1 and proceed
            v1 = v2 = 1

        # Define two dropdown widgets listing atom indices and a help box
        atom1 = widgets.Dropdown(value=v1, options=range(1, num_atoms + 1), layout=widgets.Layout(width='100px'))
        atom2 = widgets.Dropdown(value=v2, options=range(1, num_atoms + 1), layout=widgets.Layout(width='100px'))
        help_box = widgets.Button(description='?', tooltip=param['fixed_bonds']['long_glossory'], layout=widgets.Layout(width='4%'))
        box = widgets.HBox([help_box, widgets.Label(param['fixed_bonds']['glossory'], layout={'width': '400px'}), atom1, atom2])

        # Display box
        with out:
            display(box)

        # Append new widgets to the widget dictionary
        input_options['Backbone']['fixed_bonds'].append([atom1, atom2])
        

def path(file_path, param):
    """!@brief Display backbone to Jupyter notebook given a file path

    This function displays the options for the backbone. If the @a file_path does
    not exist, this function displays nothing.
    @param file_path (str) Path to a file containing the 3D structure of a backbone molecule
    @param param (dict) @a options._options_dict['Backbone']

    @returns None; @a jupyter_widgets.input_options is modified in place

    @sa upload_backbone
    @sa input_options
    @sa options._options_dict
    @sa view_nglview
    """

    if not os.path.isfile(file_path):
        # Check if the requested file is in the "pnab/data" directory
        file_path = os.path.join(__path__[0], 'data', file_path)
        if not os.path.isfile(file_path):
            # Unset the values for backbone parameters and return
            param['linker']['default'][0] = param['linker']['default'][1] = 1
            param['interconnects']['default'][0] = param['interconnects']['default'][1] = 1
            return

    # Show molecule using NGLView with atom numbers
    num_atoms = view_nglview(file_path, label=True)

    # Display widgets for backbone connection to the nucleobase and the other backbone
    linker1 = widgets.Dropdown(value=param['linker']['default'][0], options=range(1, num_atoms + 1), layout=widgets.Layout(width='100px'))
    linker2 = widgets.Dropdown(value=param['linker']['default'][1], options=range(1, num_atoms + 1), layout=widgets.Layout(width='100px'))
    interconnects1 = widgets.Dropdown(value=param['interconnects']['default'][0], options=range(1, num_atoms + 1), layout=widgets.Layout(width='100px'))
    interconnects2 = widgets.Dropdown(value=param['interconnects']['default'][1], options=range(1, num_atoms + 1), layout=widgets.Layout(width='100px'))

    # Add help boxes
    help_box = widgets.Button(description='?', tooltip=param['linker']['long_glossory'], layout=widgets.Layout(width='4%'))
    box1 = widgets.HBox([help_box, widgets.Label(param['linker']['glossory'], layout={'width': '400px'}), linker1, linker2])

    help_box = widgets.Button(description='?', tooltip=param['interconnects']['long_glossory'], layout=widgets.Layout(width='4%'))
    box2 = widgets.HBox([help_box, widgets.Label(param['interconnects']['glossory'],
                layout={'width': '400px'}), interconnects1, interconnects2])

    # Display
    display(box1)
    display(box2)

    # Add widgets to the widgets dictionary
    input_options['Backbone'] = {'file_path': file_path, 'linker': [linker1, linker2],
                           'interconnects': [interconnects1, interconnects2]}

    help_box = widgets.Button(description='?', tooltip=('This can be used to fix some of the rotatable bonds at their current values' +
                              ' to decrease the number of degrees of freedom in the search.'), layout=widgets.Layout(width='4%'))

    # Display a text widget asking for the number of fixed bonds (default to 0)
    display(widgets.HBox([help_box, widgets.interactive(fixed_bonds, num_bonds=widgets.BoundedIntText(value=len(param['fixed_bonds']['default']), min=0,
                                                        description="Number of fixed rotatable bonds", style={'description_width': 'initial'}),
                                    param=widgets.fixed(param), num_atoms=widgets.fixed(num_atoms))]))



def upload_backbone(f, param):
    """!@brief Upload a backbone file to Jupyter notebook

    This function is used to upload a backbone file and write it. If a file is
    already specified, then the backbone is displayed. If a file is uploaded, this
    function writes the file to the current working directory. Then, the function
    calls the @a jupyter_widgets.path function interactively with the backbone file path. 

    @param f (ipywidgets.FileUpload) File upload widget
    @param param (dict) @a options._options_dict['Backbone']

    @returns None

    @sa path
    @sa backbone
    """
    if f:
        # This is to unset previous values if a new file is provided
        param['linker']['default'][0] = param['linker']['default'][1] = 1
        param['interconnects']['default'][0] = param['interconnects']['default'][1] = 1

        # Get the content of the file in binary format and write it to desk
        input_file = list(f.keys())[0]
        with open(input_file, 'wb') as w:
            w.write(list(f.values())[0]['content'])

        param['file_path']['default'] = input_file

    # Add help box
    help_box = widgets.Button(description='?', tooltip=param['file_path']['long_glossory'], layout=widgets.Layout(width='4%'))

    # Use the backbone file to display the backbone options
    display(widgets.HBox([help_box, widgets.interactive(path, file_path=widgets.Text(value=param['file_path']['default'],
                                        description="Backbone File", style={'description_width': 'initial'}),
                                        param=widgets.fixed(param), upload=widgets.fixed(True))]))


def backbone(param):
    """!@brief Main backbone widget for use in Jupyter notebook

    This function displays the backbone widgets.

    @param param (dict) @a options._options_dict['Backbone']

    @returns None

    @sa upload_backbone
    @sa options._options_dict
    @sa display_options_widgets
    """

    display(widgets.HTML(value='<H3>Backbone</H3>'))

    # Display a widget for uploading backbone files
    w = widgets.interactive(upload_backbone, param=widgets.fixed(param),
                            f=widgets.FileUpload(accept='', multiple=False, 
                            description="Backbone File", style={'description_width': 'initial'}))

    # Add a help box
    help_box = widgets.Button(description='?', tooltip='Upload a backbone file here. The file will be copied to the current folder.',
                              layout=widgets.Layout(width='3%'))
    display(widgets.HBox([help_box, w]))


def base_path(file_path, param, base_number):
    """!@brief Display base to Jupyter notebook given a file path

    This function displays the options for the base. If the @a file_path does
    not exist, this function displays nothing.
    @param file_path (str) Path to a file containing the 3D structure of the base molecule
    @param param (dict) @a options._options_dict['Base']
    @param base_number The number of the new base

    @returns None; @a jupyter_widgets.input_options is modified in place

    @sa upload_base
    @sa input_options
    @sa options._options_dict
    @sa view_nglview
    """

    if not os.path.isfile(file_path):
        # Check if the requested file is in the "pnab/data" directory
        file_path = os.path.join(__path__[0], 'data', file_path)
        if not os.path.isfile(file_path):
            # Unset the values for backbone parameters and return
            param['linker']['default'][0] = param['linker']['default'][1] = 1
            return

    # Show molecule using NGLView with atom numbers
    num_atoms = view_nglview(file_path, label=True)

    # Display widgets for base-to-backbone connection
    linker1 = widgets.Dropdown(value=param['linker']['default'][0], options=range(1, num_atoms + 1), layout=widgets.Layout(width='100px'))
    linker2 = widgets.Dropdown(value=param['linker']['default'][1], options=range(1, num_atoms + 1), layout=widgets.Layout(width='100px'))

    # Add help box
    help_box = widgets.Button(description='?', tooltip=param['linker']['long_glossory'], layout=widgets.Layout(width='4%'))
    box = widgets.HBox([help_box, widgets.Label(param['linker']['glossory'], layout={'width': '400px'}), linker1, linker2])

    # Display
    display(box)

    # Code
    code = widgets.Text(value=param['code']['default'], description=param['code']['glossory'], style={'description_width': 'initial'}, layout={'width': '75%'})
    help_box = widgets.Button(description='?', tooltip=param['code']['long_glossory'], layout=widgets.Layout(width='4%'))
    box = widgets.HBox([help_box, code])
    display(box)

    # Name
    name = widgets.Text(value=param['name']['default'], description=param['name']['glossory'], style={'description_width': 'initial'}, layout={'width': '75%'})
    help_box = widgets.Button(description='?', tooltip=param['name']['long_glossory'], layout=widgets.Layout(width='4%'))
    box = widgets.HBox([help_box, name])
    display(box)

    # Pair name
    pair_name = widgets.Text(value=param['pair_name']['default'], description=param['pair_name']['glossory'], style={'description_width': 'initial'}, layout={'width': '75%'})
    help_box = widgets.Button(description='?', tooltip=param['pair_name']['long_glossory'], layout=widgets.Layout(width='4%'))
    box = widgets.HBox([help_box, pair_name])
    display(box)

    # Add widgets to the widgets dictionary
    input_options['Base %i' %base_number] = {'file_path': file_path, 'linker': [linker1, linker2], 'code': code, 'name': name, 'pair_name': pair_name}



def upload_base(f, param, base_number):
    """!@brief Upload a base file to Jupyter notebook

    This function is used to upload a base file and write it. If a file is
    already specified, then the base is displayed. If a file is uploaded, this
    function writes the file to the current working directory. Then, the function
    calls the @a jupyter_widgets.base_path function interactively with the base file path. 

    @param f (ipywidgets.FileUpload) File upload widget
    @param param (dict) @a options._options_dict['Base'] for the appropriate base number
    @param base_number (int) The number of the additional base

    @returns None

    @sa base_path
    @sa add_base
    """
    if f:
        # This is to unset previous values if a new file is provided
        param['linker']['default'][0] = param['linker']['default'][1] = 1

        # Get the content of the file in binary format and write it to desk
        input_file = list(f.keys())[0]
        with open(input_file, 'wb') as w:
            w.write(list(f.values())[0]['content'])

        param['file_path']['default'] = input_file

    # Add help box
    help_box = widgets.Button(description='?', tooltip=param['file_path']['long_glossory'], layout=widgets.Layout(width='4%'))

    # Use the base file to display the backbone options
    display(widgets.HBox([help_box, widgets.interactive(base_path, file_path=widgets.Text(value=param['file_path']['default'],
                                         description="Base File", style={'description_width': 'initial'}),
                                         param=widgets.fixed(param), upload=widgets.fixed(True), base_number=widgets.fixed(base_number))]))


def add_base(number_of_bases, param):
    """!@brief Display widgets to upload the requested number of bases

    @param number_of_bases (int) The number of additional bases to define
    @param param (dict) @a options._options_dict['Base'] and dictionaries for the other defined bases

    @sa bases
    @sa upload_base
    @sa options._options_dict
    @sa view_nglview
    """
    # Capture output. This is used to delete widget when no additional bases are requested
    out = widgets.Output()
    display(out)
    if number_of_bases == 0:
        out.clear_output()

    for i in range(number_of_bases):
        if "Base %i" %(i+1) not in param:
            param["Base %i" %(i+1)] = copy.deepcopy(param['Base'])
        # Display a widget for uploading base files
        w = widgets.interactive(upload_base, param=widgets.fixed(param["Base %i" %(i+1)]), base_number=widgets.fixed(i + 1),
                                f=widgets.FileUpload(accept='', multiple=False, 
                                description="Base File", style={'description_width': 'initial'}))

        # Add a help box
        help_box = widgets.Button(description='?', tooltip='Upload base file here. The file will be copied to the current folder.',
                                  layout=widgets.Layout(width='3%'))
        display(widgets.HBox([help_box, w]))



def bases(param):
    """!@brief Bases widget for use in Jupyter notebook

    This function displays information on the available nucleobases that are defined in the program.
    The bases are defined in "data/bases_library.yaml"

    @returns None

    @sa display_options_widgets
    """

    display(widgets.HTML(value='<H3>Bases</H3>'))
    display(widgets.HTML(value=('These bases are already defined:' +
                               '<br> Adenine (A), Guanine (G), Cytosine (C), Uracil (U), Thymine (T), ' +
                               'Cyanuric Acid (X), and Triaminopyrimidine (Y)')))
    display(widgets.HTML(value='Additional bases can also be defined. Make sure the base coordinates are in the correct reference frame.'))

    num_defined_bases = len([k for k in param if 'Base' in k]) - 1 # Subtract the item from the main options
    w = widgets.interactive(add_base, param=widgets.fixed(param), number_of_bases=widgets.BoundedIntText(value=num_defined_bases, min=0,
                                                        description="Number of additional bases", style={'description_width': 'initial'}))

    help_box = widgets.Button(description='?', tooltip=('The number of additional bases to define'), layout=widgets.Layout(width='4%'))
    display(widgets.HBox([help_box, w]))


def helical_parameters(param):
    """
    Display widgets for specifying helical parameters. It also displays an image illustrating
    the used helical parameters.

    @param param (dict) @a options._options_dict['HelicalParameters']

    @returns None; @a jupyter_widgets.input_options is modified in place

    @sa input_options
    @sa display_options_widgets 
    @sa options._options_dict
    """

    display(widgets.HTML(value='<H3>Helical Parameters</H3>'))

    # Display an image showing the helical parameters for DNA and RNA
    display(Image(os.path.join(__path__[0], 'images', 'helical_parameters.jpeg'), retina=True))

    param_dict = {}
    # Add widget to the widgets dictionary
    input_options['HelicalParameters'] = {}
    for k in param:
        param_dict[k] = []
        default = [param[k]['default'][0], param[k]['default'][1], param[k]['default'][2]] # [beginning point, end point, number of steps]
        # Set angles
        if k in ['inclination', 'tip', 'h_twist']:
            param_dict[k].append(widgets.FloatRangeSlider(value=[default[0], default[1]], min=-180, max=180, step=0.01, readout_format='.2f'))
        # Set distances
        else: # h_rise, x-displacement, y-displacement
            # limit maxium and minimum distance values to be between -10 and 10 and 0.01 Angstrom step size
            param_dict[k].append(widgets.FloatRangeSlider(value=[default[0], default[1]], min=-10, max=10, step=0.01, readout_format='.3f'))

        # Add the number of steps widgets
        param_dict[k].append(widgets.BoundedIntText(value=default[2], min=1, max=1000, step=1, description='Steps'))
        help_box = widgets.Button(description='?', tooltip=param[k]['long_glossory'], layout=widgets.Layout(width='3%'))
        box = widgets.HBox([help_box, widgets.Label(param[k]['glossory'], layout={'width': '200px'}), param_dict[k][0], param_dict[k][1]])
        display(box)

        # Add widgets for helical parameter k to the widgets dictionary
        input_options['HelicalParameters'][k] = [param_dict[k][0], param_dict[k][1]]


def algorithm(chosen_algorithm, param):
    """!@brief Display search parameters based on the chosen algorithm

    There are six search algorithms and each one has a set of input parameters.
    The search algorithms are:

    - Systematic Search: Requires specifying the dihedral angle step size, @a dihedral_step.
    - Monte Carlo Search: Requires specifying the number of steps, @a num_steps, 
        and the Monte Carlo temperature, @a monte_carlo_temperature.
    - Weighted Monte Carlo Search: Also requires specifying the weighting temperature, @a weighting_temperature.
    - Random Search: Requires specifying the number of steps, @a num_steps.
    - Weighted Random Search: Also requires specifying the weighting temperature, @a weighting_temperature.
    - Genetic Algorithm Search: Requires specifying the number of generations, @a num_steps,
        the population size, @a population_size, the mutation_rate, @a mutation_rate,
        and the crossover rate, @a crossover_rate.

    @param chosen_algorithm (str) The chosen algorithm
    @param param (dict) @a options._options_dict['RuntimeParameters']

    @returns None; @a jupyter_widgets.input_options is modified in place

    @sa runtime_parameters
    @sa input_options
    @sa options._options_dict
    """
    
    chosen_algorithm = chosen_algorithm.lower()

    # Add widget to the widgets dictionary
    input_options['RuntimeParameters']['search_algorithm'] = chosen_algorithm

    # Systematic search algorithm needs a dihedral step size
    if chosen_algorithm == 'systematic search':
        dihedral_step = widgets.BoundedFloatText(value=param['dihedral_step']['default'],
                                                 min=0.0,
                                                 max=360.0,
                                                 description=param['dihedral_step']['glossory'],
                                                 style={'description_width': 'initial'},
                                                 layout={'width': '75%'})
        help_box = widgets.Button(description='?', tooltip=param['dihedral_step']['long_glossory'], layout=widgets.Layout(width='3%'))
        display(widgets.HBox([help_box, dihedral_step]))
        input_options['RuntimeParameters']['dihedral_step'] = dihedral_step

    # The other five algorithms require specifying the number of steps or generation
    else:
        num_steps = widgets.BoundedIntText(value=param['num_steps']['default'], min=1, max=1e100,
                                           description=param['num_steps']['glossory'],
                                           style={'description_width': 'initial'},
                                           layout={'width': '75%'})
        help_box = widgets.Button(description='?', tooltip=param['num_steps']['long_glossory'], layout=widgets.Layout(width='3%'))
        display(widgets.HBox([help_box, num_steps]))
        input_options['RuntimeParameters']['num_steps'] = num_steps

    # If the algorithm weights the dihedral angle probability profile, we need to specify the temperature
    # used in the weighting
    if "weighted" in chosen_algorithm:
        weighting_temperature = widgets.BoundedFloatText(value=param['weighting_temperature']['default'], min=1, max=1e100,
                                                         description=param['weighting_temperature']['glossory'],
                                                         style={'description_width': 'initial'},
                                                         layout={'width': '75%'})
        help_box = widgets.Button(description='?', tooltip=param['weighting_temperature']['long_glossory'], layout=widgets.Layout(width='3%'))
        display(widgets.HBox([help_box, weighting_temperature]))
        input_options['RuntimeParameters']['weighting_temperature'] = weighting_temperature

    # Specify the temperature for the Monte Carlo procedure
    if "monte carlo search" in chosen_algorithm:
        monte_carlo_temperature = widgets.BoundedFloatText(value=param['monte_carlo_temperature']['default'], min=1, max=1e100,
                                                           description=param['monte_carlo_temperature']['glossory'],
                                                           style={'description_width': 'initial'},
                                                           layout={'width': '75%'})
        help_box = widgets.Button(description='?', tooltip=param['monte_carlo_temperature']['long_glossory'], layout=widgets.Layout(width='3%'))
        display(widgets.HBox([help_box, monte_carlo_temperature]))
        input_options['RuntimeParameters']['monte_carlo_temperature'] = monte_carlo_temperature

    # Specify population size, mutation rate, and crossover rate for genetic algorithm
    if 'genetic algorithm search' in chosen_algorithm:
        population_size = widgets.BoundedIntText(value=param['population_size']['default'], min=1, max=1e100,
                                                 description=param['population_size']['glossory'],
                                                 style={'description_width': 'initial'},
                                                 layout={'width': '75%'})
        help_box = widgets.Button(description='?', tooltip=param['population_size']['long_glossory'], layout=widgets.Layout(width='3%'))
        display(widgets.HBox([help_box, population_size]))
        input_options['RuntimeParameters']['population_size'] = population_size

        mutation_rate = widgets.BoundedFloatText(value=param['mutation_rate']['default'], min=0, max=1,
                                                 description=param['mutation_rate']['glossory'],
                                                 style={'description_width': 'initial'},
                                                 layout={'width': '75%'})
        help_box = widgets.Button(description='?', tooltip=param['mutation_rate']['long_glossory'], layout=widgets.Layout(width='3%'))
        display(widgets.HBox([help_box, mutation_rate]))
        input_options['RuntimeParameters']['mutation_rate'] = mutation_rate

        crossover_rate = widgets.BoundedFloatText(value=param['crossover_rate']['default'], min=0, max=1,
                                                  description=param['crossover_rate']['glossory'],
                                                  style={'description_width': 'initial'},
                                                  layout={'width': '75%'})
        help_box = widgets.Button(description='?', tooltip=param['crossover_rate']['long_glossory'], layout=widgets.Layout(width='3%'))
        display(widgets.HBox([help_box, crossover_rate]))
        input_options['RuntimeParameters']['crossover_rate'] = crossover_rate


def runtime_parameters(param):
    """!@brief Runtime parameter widget for use in Jupyter notebook

    Displays widgets for specifying runtime parameters.

    @param param (dict) @a options._options_dict['RuntimeParameters']

    @returns None; @a jupyter_widgets.input_options is modified in place

    @sa algorithm
    @sa input_options
    @sa display_options_widgets 
    @sa options._options_dict
    """

    display(widgets.HTML(value='<H3>Runtime Parameters</H3>'))

    input_options['RuntimeParameters'] = {}

    # Random number generator seed
    seed = widgets.BoundedIntText(value=param['seed']['default'], min=0, max=2**32-1,
                                  description=param['seed']['glossory'],
                                  style={'description_width': 'initial'},
                                  layout={'width': '75%'})
    help_box = widgets.Button(description='?', tooltip=param['seed']['long_glossory'], layout=widgets.Layout(width='3%'))
    display(widgets.HBox([help_box, seed]))
    input_options['RuntimeParameters']['seed'] = seed


    # Search algorithm
    display(widgets.HTML(value='<H4>Search Algorithm</H4>'))
    dropdown = widgets.Dropdown(value=param['search_algorithm']['default'].title(),
                                options=['Weighted Monte Carlo Search', 'Monte Carlo Search', 'Weighted Random Search', 'Random Search', 'Genetic Algorithm Search', 'Systematic Search'],
                                description=param['search_algorithm']['glossory'],
                                style={'description_width': 'initial'},
                                layout={'width': '75%'})
    search_algorithm = widgets.interactive_output(algorithm, {'chosen_algorithm':dropdown, 'param':widgets.fixed(param)})
    help_box = widgets.Button(description='?', tooltip=param['search_algorithm']['long_glossory'], layout=widgets.Layout(width='3%'))
    display(widgets.HBox([help_box, dropdown]))
    display(search_algorithm)

    display(widgets.HTML(value='<H4>Distance and Energy Thresholds</H4>'))

    # Distance and energy thresholds
    max_distance = widgets.BoundedFloatText(value=param['max_distance']['default'], min=0.0,
                                            description=param['max_distance']['glossory'],
                                            style={'description_width': 'initial'},
                                            layout={'width': '75%'})
    help_box = widgets.Button(description='?', tooltip=param['max_distance']['long_glossory'], layout=widgets.Layout(width='3%'))
    display(widgets.HBox([help_box, max_distance]))
    input_options['RuntimeParameters']['max_distance'] = max_distance

    # Force field
    ff_type = widgets.Dropdown(options=['GAFF', 'MMFF94', 'MMFF94s', 'UFF', 'GHEMICAL'],
                               description=param['ff_type']['glossory'],
                               style={'description_width': 'initial'},
                               layout={'width': '75%'})
    help_box = widgets.Button(description='?', tooltip=param['ff_type']['long_glossory'], layout=widgets.Layout(width='3%'))
    display(widgets.HBox([help_box, ff_type]))
    input_options['RuntimeParameters']['ff_type'] = ff_type

    # Energy filters
    input_options['RuntimeParameters']['energy_filter'] = []
    for i in range(5):
        label = param['energy_filter']['glossory'].split('\n')[i]
        energy_filter = widgets.FloatText(value=param['energy_filter']['default'][i],
                                          description=label,
                                          style={'description_width': 'initial'},
                                          layout={'width': '75%'})
        help_box = widgets.Button(description='?', tooltip=param['energy_filter']['long_glossory'].split('\n')[i], layout=widgets.Layout(width='3%'))
        display(widgets.HBox([help_box, energy_filter]))
        input_options['RuntimeParameters']['energy_filter'].append(energy_filter)


    display(widgets.HTML(value='<H4>Structural Parameters</H4>'))
    # Base sequence
    strand = widgets.Text(value=''.join(param['strand']['default']),
                          description=param['strand']['glossory'],
                          style={'description_width': 'initial'},
                          layout={'width': '75%'})
    help_box = widgets.Button(description='?', tooltip=param['strand']['long_glossory'], layout=widgets.Layout(width='3%'))
    display(widgets.HBox([help_box, strand]))
    input_options['RuntimeParameters']['strand'] = strand

    # Pair adenine with uracil? Default is A-T base pair
    pair_A_U = widgets.Checkbox(value=param['pair_A_U']['default'], indent=False,
                                description=param['pair_A_U']['glossory'],
                                style={'description_width': 'initial'},
                                layout={'width': '75%'})
    help_box = widgets.Button(description='?', tooltip=param['pair_A_U']['long_glossory'], layout=widgets.Layout(width='3%'))
    display(widgets.HBox([help_box, pair_A_U]))
    input_options['RuntimeParameters']['pair_A_U'] =  pair_A_U

    # Is hexad
    is_hexad = widgets.Checkbox(value=param['is_hexad']['default'], indent=False,
                                description=param['is_hexad']['glossory'],
                                style={'description_width': 'initial'},
                                layout={'width': '75%'})
    help_box = widgets.Button(description='?', tooltip=param['is_hexad']['long_glossory'], layout=widgets.Layout(width='3%')) 
    display(widgets.HBox([help_box, is_hexad]))
    input_options['RuntimeParameters']['is_hexad'] =  is_hexad

    # Build strand
    input_options['RuntimeParameters']['build_strand'] = []
    help_box = widgets.Button(description='?', tooltip=param['build_strand']['long_glossory'], layout=widgets.Layout(width='3%'))
    box = [help_box, widgets.Label(param['build_strand']['glossory'])]
    for i in range(6):
        build_strand = widgets.Checkbox(value=param['build_strand']['default'][i], indent=False, layout={'width': '50px'})
        input_options['RuntimeParameters']['build_strand'].append(build_strand)
        box.append(build_strand)
    box = widgets.HBox(box, layout={'width':'100%'})
    display(box)

    # Orientation of each strand in the hexad
    input_options['RuntimeParameters']['strand_orientation'] = []
    help_box = widgets.Button(description='?', tooltip=param['strand_orientation']['long_glossory'], layout=widgets.Layout(width='3%'))
    box = [help_box, widgets.Label(param['strand_orientation']['glossory'])]
    for i in range(6):
        strand_orientation = widgets.Checkbox(value=param['strand_orientation']['default'][i], indent=False, layout={'width': '50px'})
        input_options['RuntimeParameters']['strand_orientation'].append(strand_orientation)
        box.append(strand_orientation)
    box = widgets.HBox(box, layout={'width':'100%'})
    display(box)


def upload_input(param, f):
    """!@brief Widget to upload an input file

    If a file is uploaded, it is written in the current working directory.
    Then, this function calls @a jupyter_widgets.display_options_widgets function with the newly
    written file.

    @param param (dict) @a options._options_dict
    @param f (ipywidgets.FileUpload) file upload widget

    @returns None

    @sa display_options_widgets
    """
    if f:
        # If an input file is provided, write it to the current working directory
        input_file = list(f.keys())[0]
        with open(input_file, 'wb') as w:
            w.write(list(f.values())[0]['content'])

        # Display the name of the input file
        display(widgets.HTML(value=input_file))

        # Display all the widgets with the options specified in the input file
        display_options_widgets(param, input_file, uploaded=True)


def display_options_widgets(param, input_file, uploaded=False):
    """!@brief Display all widgets in Jupyter notebook given an input file

    If @a input_file is provided, then it updates the default options in
    @a options._options_dict to display the user-defined options. If @a
    input_file is "Upload file", then an ipywidgets.FileUpload widget is displayed.
    This function also displays a widget for running the program after the
    user defines all the options.


    @param param (dict) @a options._options_dict
    @param input_file (str) Input file
    @param uploaded (bool) Whether a file is uploaded by the user or not

    @returns None

    @sa backbone
    @sa bases
    @sa helical_parameters
    @sa runtime_parameters
    @sa upload_input
    @sa run
    @sa user_input_file
    @sa options._options_dict
    """

    if input_file == "Upload file":
        # Display a widget for uploading input file
        w = widgets.interactive(upload_input, param=widgets.fixed(param), f=widgets.FileUpload(accept='', multiple=False, description="Input File"))
        help_box = widgets.Button(description='?', tooltip='Upload your input file here. The file will be copied to the current folder.', layout=widgets.Layout(width='4%'))
        display(widgets.HBox([help_box, w]))
        return

    # If the file is not uploaded, then we are using one of the example files
    # in the pnab/data directory
    if not uploaded:
        input_file = os.path.join(__path__[0], 'data', input_file)

    options = yaml.load(open(input_file, 'r'), yaml.FullLoader)

    # Clean dictionary from additional bases if they are not defined in the options
    # This is necessary when the user switches from an input file that has additional bases to one that does not
    num_defined_bases = len([k for k in param if 'Base' in k]) - 1 # Subtract the item from the main options
    for i in range(1, num_defined_bases + 1):
        param.pop('Base %i' %i)
        input_options.pop('Base %i' %i)

    # Update the default options to display those provided by the user
    for k1, v1 in options.items():
        if 'Base' in k1:
            param[k1] = copy.deepcopy(param['Base'])
        for k2, v2 in options[k1].items():
            param[k1][k2]['default'] = options[k1][k2]

    # Display all options widgets
    backbone(param['Backbone'])
    bases_param = {k:val for k, val in param.items() if 'Base' in k}
    bases(bases_param)
    helical_parameters(param['HelicalParameters'])
    runtime_parameters(param['RuntimeParameters'])

    # Display run widget
    button = widgets.Button(description='Run', tooltip='Click here to run the program with the provided options. Once the program finishes, the results will be displayed below.')
    button.on_click(run)
    display(button)


def user_input_file(param):
    """!@brief Display input file options

    Gives the user several input files that they can try, and the option that
    they upload their input file. Once an input file is chosen, all widgets
    are displayed with the user-defined options.

    @param param (dict) @a options._options_dict

    @returns None

    @sa display_options_widgets
    """

    # Provide three input files as examples
    w = widgets.interactive(display_options_widgets, param=widgets.fixed(param), uploaded=widgets.fixed(False),
            input_file=widgets.Dropdown(options=['RNA.yaml', 'DNA.yaml', 'Hexad.yaml', 'Upload file'],
                                        style={'description_width': 'initial'}, description='Input File'))
    help_box = widgets.Button(description='?', tooltip=('There are existing example files for RNA, DNA, and hexad geometries.' +
                                                        ' You can use these examples as a starting point for customizing your input options.' + 
                                                        ' Alternatively, you can upload your own input file.'),
               layout=widgets.Layout(width='3%'))
    display(widgets.HBox([help_box, w]))


def run(button):
    """!@brief Function to run the code when the user finishes specifying all options

    It extracts the user-defined options and creates a @a pNAB.pNAB instance. Then,
    it run the code and display the results.

    @param button (variable) variable for using the widgets.Button.on_click method

    @sa display_options_widgets
    @sa pNAB.pNAB
    @sa show_results
    """

    # Extract options
    run_options = extract_options()

    # Write a file named "options.yaml" containing all the user-defined options
    file_path = 'options.yaml'

    # Make sure we do not overwrite other files that have the same name
    # by prepending enough "_"
    while True:
        if os.path.isfile(file_path):
            file_path = '_' + file_path
        else:
            if os.path.isfile('options.yaml'):
                os.rename('options.yaml', file_path)
            break

    # Write time stamps
    time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    with open('options.yaml', 'w') as f:
        f.write('# ' + time + '\n')
        f.write(yaml.dump(run_options))

    # Capture progress report 
    out = widgets.Output()
    display(out)

    with out:
        # Run the code
        run = pNAB('options.yaml')
        run.run()

    # Delete progress report
    out.clear_output()

    # Get output files
    files = [str(int(conformer[0])) + '_' + str(int(conformer[1])) + '.pdb' for conformer in run.results]

    time = datetime.datetime.now().strftime('%Y%m%d%H%M%S')
    with ZipFile('output' + time + '.zip', 'w') as z:
        for f in ['options.yaml', 'results.csv', 'summary.csv', 'prefix.yaml']:
            z.write(f)
        for f in files:
            z.write(f)

    display(widgets.HTML("""<a href="output""" + time + """.zip" target="_blank">Download Output</a>"""))
    display(Javascript("""var url="output%s.zip"\nwindow.open(url, 'download')""" %time))

    # If no results are found, print and return
    if run.results.size == 0:
        print("No candidate found")
    # else display conformers and their properties
    else:
        # Sort by total energy
        run.results = run.results[run.results[:, 7].argsort()]
        show_results(run.results, run.header, run.prefix)



def extract_options():
    """!@brief Extracts user options from the widgets

    This function extracts the values of all options specifed by the 
    user. It access all the widgets in @a jupyter_widgets.input_options
    and extracts their values.

    @returns user_options A dictionary of the values of all the user-defined options

    @sa run
    @sa input_options
    @sa options._options_dict
    """

    user_options = {}
    # I did not use a smart way to get the values, but this should work
    for k1, val1 in input_options.items():
        user_options[k1] = {}
        for k2, val2 in val1.items():
            if isinstance(val2, str):
                user_options[k1][k2] = val2
            elif k1 == 'Backbone' or 'Base' in k1:
                if k2 == 'fixed_bonds':
                    user_options[k1][k2] = []
                    for l in val2:
                        user_options[k1][k2].append([l[0].value, l[1].value])
                elif isinstance(val2, list):
                    user_options[k1][k2] = [val2[0].value, val2[1].value]
                elif isinstance(val2, str):
                    user_options[k1][k2] = val2
                else:
                    user_options[k1][k2] = val2.value

            elif k1 == 'HelicalParameters':
                user_options[k1][k2] = [val2[0].value[0], val2[0].value[1], val2[1].value]

            else:
                if k2 == 'energy_filter' or k2 == 'strand_orientation' or k2 == 'build_strand':
                    user_options[k1][k2] = [i.value for i in val2]

                else:
                    user_options[k1][k2] = val2.value

    return user_options



def single_result(result, header, results, prefix):
    """!@brief Interactive function to display a single result

    It displays one accepted conformer and prints its energies and other properties.

    @param result (int) index of the conformer in the @a results array
    @param header (str) A comma separated string of the output properties of the conformer
    @param results (np.ndarray) a numpy array of all the accepted conformers
    @param prefix (dict) A dictionary of the prefix of the run and the associated helical configuration

    @returns None    

    @sa show_results
    """

    # Extract result from the array
    result = results[result]

    # Get the name of the PDB file containing the conformer
    conformer = str(int(result[0])) + '_' + str(int(result[1])) + '.pdb'

    # Print information
    print(conformer)
    print(prefix['%i' %result[0]])

    # Display a color based on energy terms
    cost = sum([result[i] for i in range(3, 7)])
    color = 'green' if cost < 2 else 'yellow' if cost < 5 else 'red'
    display(widgets.ColorPicker(concise=True, value=color, disabled=True))

    for i in range(2, len(result)):
        print(header.split(', ')[i] + ': %.3f' %result[i])

    # Display the molecule
    view_nglview(conformer)
    

def show_results(results, header, prefix):
    """!@brief Display results

    A function to display all the accepted candidates using a dropdown list.

    @param results (np.ndarray) numpy array of results sorted by total energy
    @param header (str) A comma separated string of the output properties of the conformer
    @param prefix (dict) A dictionary of the prefix of the run and the associated helical configuration

    @returns None

    @sa run
    @sa single_result 
    """

    # Set a list of conformer names and their indices in the results
    options = [(str(int(conformer[0])) + '_' + str(int(conformer[1])) + '.pdb', i) for i, conformer in enumerate(results)]

    # Show a dropdown widget of all the accepted conformers
    dropdown = widgets.Dropdown(value=options[0][1], options=options, style={'description_width': 'initial'}, description='Conformer')
    w = widgets.interactive(single_result, result=dropdown, header=widgets.fixed(header), results=widgets.fixed(results), prefix=widgets.fixed(prefix))
    display(w)



def builder():
    """!@brief The function called from the Jupyter notebook to display the widgets

    Execute the following to display all the widgets in the notebook
    @code
    import pnab
    pnab.builder()
    @endcode

    @returns None
    """

    # Prevents auto-scrolling of the notebook
    disable_js = """
IPython.OutputArea.prototype._should_scroll = function(lines) {
    return false;
}
"""
    display(Javascript(disable_js))

    user_input_file(_options_dict)
