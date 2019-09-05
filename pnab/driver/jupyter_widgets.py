import os
import datetime

import yaml
import ipywidgets as widgets
from IPython.display import display, Javascript

from pnab.driver.pNAB import pNAB
from pnab.driver import draw
from pnab.driver.options import _options_dict as _options
from pnab import __path__


options = {}

def fixed_bonds(n, num_atoms, param):
    options['Backbone']['fixed_bonds'] = []
    out = widgets.Output()
    display(out)
    if n == 0:
        out.clear_output()
    for i in range(n):
        try:
            v1, v2 = param['fixed_bonds']['default'][i][0], param['fixed_bonds']['default'][i][1]
        except IndexError:
            v1 = v2 = 1
        atom1 = widgets.Dropdown(value=v1, options=range(1, num_atoms + 1), layout=widgets.Layout(width='100px'))
        atom2 = widgets.Dropdown(value=v2, options=range(1, num_atoms + 1), layout=widgets.Layout(width='100px'))
        help_box = widgets.Button(description='?', tooltip=param['fixed_bonds']['long_glossory'], layout=widgets.Layout(width='4%'))
        box = widgets.HBox([help_box, widgets.Label(param['fixed_bonds']['glossory'], layout={'width': '400px'}), atom1, atom2])
        with out:
            display(box)
        options['Backbone']['fixed_bonds'].append([atom1, atom2])
        

def path(file_path, param):
    """Display backbone to Jupyter notebook given a file path"""

    if not os.path.isfile(file_path):
        # Check if the requested file is in the "pnab/data" directory
        file_path = os.path.join(__path__[0], 'data', file_path)
        if not os.path.isfile(file_path):
            # Unset the values for backbone parameters and return
            param['linker']['default'][0] = param['linker']['default'][1] = 1
            param['interconnects']['default'][0] = param['interconnects']['default'][1] = 1
            return

    # Show molecule using NGLView with atom numbers
    num_atoms = draw.view_nglview(file_path, label=True)

    # Display widgets for backbone connection to the nucleobase and the other backbone
    linker1 = widgets.Dropdown(value=param['linker']['default'][0], options=range(1, num_atoms + 1), layout=widgets.Layout(width='100px'))
    linker2 = widgets.Dropdown(value=param['linker']['default'][1], options=range(1, num_atoms + 1), layout=widgets.Layout(width='100px'))
    interconnects1 = widgets.Dropdown(value=param['interconnects']['default'][0], options=range(1, num_atoms + 1), layout=widgets.Layout(width='100px'))
    interconnects2 = widgets.Dropdown(value=param['interconnects']['default'][1], options=range(1, num_atoms + 1), layout=widgets.Layout(width='100px'))

    help_box = widgets.Button(description='?', tooltip=param['linker']['long_glossory'], layout=widgets.Layout(width='4%'))
    box1 = widgets.HBox([help_box, widgets.Label(param['linker']['glossory'], layout={'width': '400px'}), linker1, linker2])

    help_box = widgets.Button(description='?', tooltip=param['interconnects']['long_glossory'], layout=widgets.Layout(width='4%'))
    box2 = widgets.HBox([help_box, widgets.Label(param['interconnects']['glossory'],
                layout={'width': '400px'}), interconnects1, interconnects2])

    display(box1)
    display(box2)

    options['Backbone'] = {'file_path': file_path, 'linker': [linker1, linker2],
                           'interconnects': [interconnects1, interconnects2]}

    help_box = widgets.Button(description='?', tooltip=('This can be used to fix some of the rotatable bonds at their current values' +
                              ' to decrease the number of degrees of freedom in the search.'), layout=widgets.Layout(width='4%'))
    display(widgets.HBox([help_box, widgets.interactive(fixed_bonds, n=widgets.BoundedIntText(value=len(param['fixed_bonds']['default']), min=0,
                                                        description="Number of fixed rotatable bonds", style={'description_width': 'initial'}),
                                    param=widgets.fixed(param), num_atoms=widgets.fixed(num_atoms))]))



def upload_backbone(f, param):
    """Upload a backbone file to Jupyter notebook"""
    if f:
        param['linker']['default'][0] = param['linker']['default'][1] = 1
        param['interconnects']['default'][0] = param['interconnects']['default'][1] = 1

        # Get the content of the file in binary format and write it to desk
        input_file = list(f.keys())[0]
        with open(input_file, 'wb') as w:
            w.write(list(f.values())[0]['content'])

        param['file_path']['default'] = input_file

    help_box = widgets.Button(description='?', tooltip=param['file_path']['long_glossory'], layout=widgets.Layout(width='4%'))
    display(widgets.HBox([help_box, widgets.interactive(path, file_path=widgets.Text(value=param['file_path']['default'],
                                        description="Backbone File", style={'description_width': 'initial'}),
                                        param=widgets.fixed(param), upload=widgets.fixed(True))]))


def backbone(param):
    """Backbone widget for use in Jupyter notebook"""

    display(widgets.HTML(value='<b>Backbone</b>'))

    display(widgets.HTML(value='Draw a backbone molecule if you need then click on "done!"'))
    draw.draw()

    w = widgets.interactive(upload_backbone, param=widgets.fixed(param),
                            f=widgets.FileUpload(accept='', multiple=False, 
                            description="Backbone File", style={'description_width': 'initial'}))
    help_box = widgets.Button(description='?', tooltip='Upload a backbone file here. The file will be copied to the current folder.',
                              layout=widgets.Layout(width='3%'))
    display(widgets.HBox([help_box, w]))


def bases():
    """Bases widget for use in Jupyter notebook"""

    display(widgets.HTML(value='<b>Bases</b>'))
    display(widgets.HTML(value=('These bases are already defined:' +
                               '<br> Adenine (A), Guanine (G), Cytosine (C), Uracil (U), Thymine (T), ' +
                               'Cyanuric Acid (X), and Triaminopyrimidine (Y)')))



def helical_parameters(param):
    """Helical parameter widget for use in Jupyter notebook"""

    display(widgets.HTML(value='<b>Helical Parameters</b>'))

    param_dict = {}
    options['HelicalParameters'] = {}
    for k in param:
        param_dict[k] = []
        default = [param[k]['default'][0], param[k]['default'][1], param[k]['default'][2]] # [beginning point, end point, number of steps]
        # Set angles
        if k in ['inclination', 'tip', 'h_twist']:
            param_dict[k].append(widgets.FloatRangeSlider(value=[default[0], default[1]], min=-180, max=180, step=0.01, readout_format='.2f'))
        # Set distances
        else:
            param_dict[k].append(widgets.FloatRangeSlider(value=[default[0], default[1]], min=-10, max=10, step=0.01, readout_format='.3f'))
        param_dict[k].append(widgets.BoundedIntText(value=default[2], min=1, max=1000, step=1, description='Steps'))
        help_box = widgets.Button(description='?', tooltip=param[k]['long_glossory'], layout=widgets.Layout(width='3%'))
        box = widgets.HBox([help_box, widgets.Label(param[k]['glossory'], layout={'width': '200px'}), param_dict[k][0], param_dict[k][1]])
        display(box)
        options['HelicalParameters'][k] = [param_dict[k][0], param_dict[k][1]]


def algorithm(f, param):
    """Display search parameters based on the chosen algorithm"""

    options['RuntimeParameters']['search_algorithm'] = f.lower()
    if 'random search' in f.lower() or "monte carlo search" in f.lower() or 'genetic algorithm search' in f.lower():
        num_steps = widgets.IntText(value=param['num_steps']['default'],
                                    description=param['num_steps']['glossory'],
                                    style={'description_width': 'initial'},
                                    layout={'width': '75%'})
        help_box = widgets.Button(description='?', tooltip=param['num_steps']['long_glossory'], layout=widgets.Layout(width='3%'))
        display(widgets.HBox([help_box, num_steps]))
        options['RuntimeParameters']['num_steps'] = num_steps

    elif f.lower() == 'systematic search':
        dihedral_step = widgets.BoundedFloatText(value=param['dihedral_step']['default'],
                                                 min=0.0,
                                                 max=360.0,
                                                 description=param['dihedral_step']['glossory'],
                                                 style={'description_width': 'initial'},
                                                 layout={'width': '75%'})
        help_box = widgets.Button(description='?', tooltip=param['dihedral_step']['long_glossory'], layout=widgets.Layout(width='3%'))
        display(widgets.HBox([help_box, dihedral_step]))
        options['RuntimeParameters']['dihedral_step'] = dihedral_step


def runtime_parameters(param):
    """Runtime parameter widget for use in Jupyter notebook"""

    display(widgets.HTML(value='<b>Runtime Parameters</b>'))

    options['RuntimeParameters'] = {}

    # Search algorithm
    dropdown = widgets.Dropdown(value=param['search_algorithm']['default'].title(),
                                options=['Weighted Monte Carlo Search', 'Monte Carlo Search', 'Weighted Random Search', 'Random Search', 'Genetic Algorithm Search', 'Systematic Search'],
                                description=param['search_algorithm']['glossory'],
                                style={'description_width': 'initial'},
                                layout={'width': '75%'})
    search_algorithm = widgets.interactive_output(algorithm, {'f':dropdown, 'param':widgets.fixed(param)})
    help_box = widgets.Button(description='?', tooltip=param['search_algorithm']['long_glossory'], layout=widgets.Layout(width='3%'))
    display(widgets.HBox([help_box, dropdown]))
    display(search_algorithm)

    # Force field
    ff_type = widgets.Dropdown(options=['GAFF', 'MMFF94', 'MMFF94s', 'UFF', 'GHEMICAL'],
                               description=param['type']['glossory'],
                               style={'description_width': 'initial'},
                               layout={'width': '75%'})
    help_box = widgets.Button(description='?', tooltip=param['type']['long_glossory'], layout=widgets.Layout(width='3%'))
    display(widgets.HBox([help_box, ff_type]))
    options['RuntimeParameters']['type'] = ff_type

    # Distance and energy thresholds
    max_distance = widgets.FloatText(value=param['max_distance']['default'],
                                     description=param['max_distance']['glossory'],
                                     style={'description_width': 'initial'},
                                     layout={'width': '75%'})
    help_box = widgets.Button(description='?', tooltip=param['max_distance']['long_glossory'], layout=widgets.Layout(width='3%'))
    display(widgets.HBox([help_box, max_distance]))
    options['RuntimeParameters']['max_distance'] = max_distance

    options['RuntimeParameters']['energy_filter'] = []
    for i in range(5):
        label = param['energy_filter']['glossory'].split('\n')[i]
        energy_filter = widgets.FloatText(value=param['energy_filter']['default'][i],
                                          description=label,
                                          style={'description_width': 'initial'},
                                          layout={'width': '75%'})
        help_box = widgets.Button(description='?', tooltip=param['energy_filter']['long_glossory'].split('\n')[i], layout=widgets.Layout(width='3%'))
        display(widgets.HBox([help_box, energy_filter]))
        options['RuntimeParameters']['energy_filter'].append(energy_filter)

    # Base sequence
    strand = widgets.Text(value=''.join(param['strand']['default']),
                          description=param['strand']['glossory'],
                          style={'description_width': 'initial'},
                          layout={'width': '75%'})
    help_box = widgets.Button(description='?', tooltip=param['strand']['long_glossory'], layout=widgets.Layout(width='3%'))
    display(widgets.HBox([help_box, strand]))
    options['RuntimeParameters']['strand'] = strand

    # Is double stranded
    is_double_stranded = widgets.Checkbox(value=param['is_double_stranded']['default'], indent=False,
                                          description=param['is_double_stranded']['glossory'],
                                          style={'description_width': 'initial'},
                                          layout={'width': '75%'})
    help_box = widgets.Button(description='?', tooltip=param['is_double_stranded']['long_glossory'], layout=widgets.Layout(width='3%'))
    display(widgets.HBox([help_box, is_double_stranded]))
    options['RuntimeParameters']['is_double_stranded'] =  is_double_stranded

    # Pair adenine with uracil? Default is A-T base pair
    pair_A_U = widgets.Checkbox(value=param['pair_A_U']['default'], indent=False,
                                description=param['pair_A_U']['glossory'],
                                style={'description_width': 'initial'},
                                layout={'width': '75%'})
    help_box = widgets.Button(description='?', tooltip=param['pair_A_U']['long_glossory'], layout=widgets.Layout(width='3%'))
    display(widgets.HBox([help_box, pair_A_U]))
    options['RuntimeParameters']['pair_A_U'] =  pair_A_U

    # Is hexad
    is_hexad = widgets.Checkbox(value=param['is_hexad']['default'], indent=False,
                                description=param['is_hexad']['glossory'],
                                style={'description_width': 'initial'},
                                layout={'width': '75%'})
    help_box = widgets.Button(description='?', tooltip=param['is_hexad']['long_glossory'], layout=widgets.Layout(width='3%')) 
    display(widgets.HBox([help_box, is_hexad]))
    options['RuntimeParameters']['is_hexad'] =  is_hexad

    # Orientation of each strand in the hexad
    options['RuntimeParameters']['strand_orientation'] = []
    help_box = widgets.Button(description='?', tooltip=param['strand_orientation']['long_glossory'], layout=widgets.Layout(width='3%'))
    box = [help_box, widgets.Label(param['strand_orientation']['glossory'])]
    for i in range(6):
        strand_orientation = widgets.Checkbox(value=param['strand_orientation']['default'][i], indent=False, layout={'width': '50px'})
        options['RuntimeParameters']['strand_orientation'].append(strand_orientation)
        box.append(strand_orientation)
    box = widgets.HBox(box, layout={'width':'100%'})
    display(box)


def upload_input(param, f):
    """widgets to upload an input file"""
    if f:
        input_file = list(f.keys())[0]
        with open(input_file, 'wb') as w:
            w.write(list(f.values())[0]['content'])

        display(widgets.HTML(value=input_file))
        show(param, input_file, uploaded=True)


def show(param, input_file, uploaded=False):
    """Display all widgets in Jupyter notebook"""

    if input_file == "Upload file":
        w = widgets.interactive(upload_input, param=widgets.fixed(param), f=widgets.FileUpload(accept='', multiple=False, description="Input File"))
        help_box = widgets.Button(description='?', tooltip='Upload your input file here. The file will be copied to the current folder.', layout=widgets.Layout(width='4%'))
        display(widgets.HBox([help_box, w]))
        return

    if not uploaded:
        input_file = os.path.join(__path__[0], 'data', input_file)

    options = yaml.load(open(input_file, 'r'), yaml.FullLoader)

    for k1, v1 in options.items():
        for k2, v2 in options[k1].items():
            param[k1][k2]['default'] = options[k1][k2]

    backbone(param['Backbone'])
    bases()
    helical_parameters(param['HelicalParameters'])
    runtime_parameters(param['RuntimeParameters'])


def run(b):
    """Function to run the code when the user finishes specifying all options"""

    run_options = extract_options(options)
    file_path = 'options.yaml'
    while True:
        if os.path.isfile(file_path):
            file_path = '_' + file_path
        else:
            if os.path.isfile('options.yaml'):
                os.rename('options.yaml', file_path)
            break

    time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    with open('options.yaml', 'w') as f:
        f.write('# ' + time + '\n')
        f.write(yaml.dump(run_options))

    run = pNAB('options.yaml')
    run.run()
    if run.results.size == 0:
        print("No candidate found")
    else:
        show_results(run.results, run.header, run.prefix)


def display_widgets(param):
    """Display all widgets Jupyter notebook and options

    return widgets which can be used to extract user options
    """

    # Prevents auto-scrolling of the notebook
    disable_js = """
IPython.OutputArea.prototype._should_scroll = function(lines) {
    return false;
}
"""
    display(Javascript(disable_js))
    # Provide three input files as examples
    w = widgets.interactive(show, param=widgets.fixed(param), uploaded=widgets.fixed(False),
            input_file=widgets.Dropdown(options=['RNA.yaml', 'DNA.yaml', 'Hexad.yaml', 'Upload file'],
                                        style={'description_width': 'initial'}, description='Input File'))
    help_box = widgets.Button(description='?', tooltip=('There are existing example files for RNA, DNA, and hexad geometries.' +
                                                        ' You can use these examples as a starting point for customizing your input options.' + 
                                                        ' Alternatively, you can upload your own input file.'),
               layout=widgets.Layout(width='3%'))
    display(widgets.HBox([help_box, w]))
    # Run widget
    button = widgets.Button(description='Run', tooltip='Click here to run the program with the provided options. Once the program finishes, the results will be displayed below.')
    button.on_click(run)
    display(button)


def extract_options(param):
    """Extract user options from the widgets"""

    user_options = {}
    for k1, val1 in param.items():
        user_options[k1] = {}
        for k2, val2 in val1.items():
            if isinstance(val2, str):
                user_options[k1][k2] = val2
            elif k1 == 'Backbone':
                if k2 == 'fixed_bonds':
                    user_options[k1][k2] = []
                    for l in val2:
                        user_options[k1][k2].append([l[0].value, l[1].value])
                elif isinstance(val2, list):
                    user_options[k1][k2] = [val2[0].value, val2[1].value]
                else:
                    user_options[k1][k2] = val2

            elif k1 == 'HelicalParameters':
                user_options[k1][k2] = [val2[0].value[0], val2[0].value[1], val2[1].value]

            else:
                if k2 == 'energy_filter' or k2 == 'strand_orientation':
                    user_options[k1][k2] = [i.value for i in val2]

                else:
                    user_options[k1][k2] = val2.value

    return user_options


def builder():
    """Funtion called from the Jupyter notebook to display input"""

    display_widgets(_options)

def single_result(result, header, results, prefix):
    """Interactive function to display a single result"""

    result = results[result]
    conformer = str(int(result[0])) + '_' + str(int(result[1])) + '.pdb'
    draw.view_nglview(conformer)
    print(conformer)
    print(prefix['%i' %result[0]])
    for i in range(2, len(result)):
        print(header.split(', ')[i] + ': ' + str(result[i]))
    

def show_results(results, header, prefix):
    """Display results"""

    options = [(str(int(conformer[0])) + '_' + str(int(conformer[1])) + '.pdb', i) for i, conformer in enumerate(results)]
    dropdown = widgets.Dropdown(value=options[0][1], options=options, style={'description_width': 'initial'}, description='Conformer')
    w = widgets.interactive(single_result, result=dropdown, header=widgets.fixed(header), results=widgets.fixed(results), prefix=widgets.fixed(prefix))
    display(w)
