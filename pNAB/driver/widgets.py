import ipywidgets as widgets
from IPython.display import display

from pNAB.driver import draw


options = {}

def backbone(param, file_path):
    """Backbone widget for use in Jupyter notebook"""

    import os

    if not os.path.isfile(file_path):
        return

    num_atoms = draw.view_py3dmol(file_path, label=True)

    linker = widgets.SelectMultiple(options=range(1, num_atoms + 1)) 
    interconnects = widgets.SelectMultiple(options=range(1, num_atoms + 1))

    box1 = widgets.HBox([widgets.Label(param['linker']['glossory'], layout={'width': '400px'}), linker])
    box2 = widgets.HBox([widgets.Label(param['interconnects']['glossory'], layout={'width': '400px'}), interconnects])

    display(box1)
    display(box2)

    options['Backbone'] = {'file_path': file_path, 'linker': linker,
                           'interconnects': interconnects}


def base(param, file_path, base_number):
    """Base widget for use in Jupyter notebook"""
    import os

    if not os.path.isfile(file_path):
        return

    num_atoms = draw.view_py3dmol(file_path, label=True)

    linker = widgets.SelectMultiple(options=range(1, num_atoms + 1)) 
    code = widgets.Text()
    name = widgets.Text()
    pair_name = widgets.Text()

    box1 = widgets.HBox([widgets.Label(param['linker']['glossory'], layout={'width': '400px'}), linker])
    box2 = widgets.HBox([widgets.Label(param['code']['glossory'], layout={'width': '400px'}), code])
    box3 = widgets.HBox([widgets.Label(param['name']['glossory'], layout={'width': '400px'}), name])
    box4 = widgets.HBox([widgets.Label(param['pair_name']['glossory'], layout={'width': '400px'}), pair_name])

    display(box1)
    display(box2)
    display(box3)
    display(box4)

    options['Base %i' %base_number]['file_path'] = file_path
    options['Base %i' %base_number]['linker'] = linker
    options['Base %i' %base_number]['code'] = code
    options['Base %i' %base_number]['name'] = name
    options['Base %i' %base_number]['pair_name'] = pair_name


def bases(param, num_bases):
    """Number of bases widget for use in Jupyter notebook

    Dynamically changes the number of base options based on the number of bases
    """

    for i in range(num_bases):
        options['Base %i' %(i+1)] = {}
        path = widgets.interactive(base, param=widgets.fixed(param), file_path=widgets.Text(value='', description="Base File", style={'description_width': 'initial'}),
                                   base_number=widgets.fixed(i+1))
        display(path)


def helical_parameters(param):
    """Helical parameter widget for use in Jupyter notebook"""

    param_dict = {}
    options['HelicalParameters'] = {}
    for k in param:
        param_dict[k] = []
        if k in ['inclination', 'tip', 'twist']:
            param_dict[k].append(widgets.FloatRangeSlider(value=[-180, 180], min=-180, max=180, step=0.01, readout_format='.2f'))
        else:
            param_dict[k].append(widgets.FloatRangeSlider(value=[-10, 100], min=-10, max=10, step=0.01, readout_format='.3f'))
        param_dict[k].append(widgets.BoundedIntText(value=1, min=1, max=1000, step=1, description='Steps'))
        box = widgets.HBox([widgets.Label(param[k]['glossory'], layout={'width': '100px'}), param_dict[k][0], param_dict[k][1]])
        display(box)
        options['HelicalParameters'][k] = [param_dict[k][0], param_dict[k][1]]


def runtime_parameters(param):
    """Runtime parameter widget for use in Jupyter notebook"""

    options['RuntimeParameters'] = {}

    num_steps = widgets.IntText(value=param['num_steps']['default'])
    box = widgets.HBox([widgets.Label(param['num_steps']['glossory'], layout={'width': '400px'}), num_steps])
    display(box)
    options['RuntimeParameters']['num_steps'] = num_steps

    ff_type = widgets.Text(value=param['type']['default'])
    box = widgets.HBox([widgets.Label(param['type']['glossory'], layout={'width': '400px'}), ff_type])
    display(box)
    options['RuntimeParameters']['type'] = ff_type

    algorithm = widgets.Text(value=param['algorithm']['default'])
    box = widgets.HBox([widgets.Label(param['algorithm']['glossory'], layout={'width': '400px'}), algorithm])
    display(box)
    options['RuntimeParameters']['algorithm'] = algorithm

    options['RuntimeParameters']['energy_filter'] = []
    for i in range(5): 
        label = param['energy_filter']['glossory'].split('\n')[i]
        energy_filter = widgets.FloatText(value=param['energy_filter']['default'][i])
        box = widgets.HBox([widgets.Label(label, layout={'width': '400px'}), energy_filter])
        display(box)
        options['RuntimeParameters']['energy_filter'].append(energy_filter)

    max_distance = widgets.FloatText(value=param['max_distance']['default'])
    box = widgets.HBox([widgets.Label(param['max_distance']['glossory'], layout={'width': '400px'}), max_distance])
    display(box)
    options['RuntimeParameters']['max_distance'] = max_distance

    strand = widgets.Text(value=param['strand']['default'])
    box = widgets.HBox([widgets.Label(param['strand']['glossory'], layout={'width': '400px'}), strand])
    display(box)
    options['RuntimeParameters']['strand'] = strand

    is_double_stranded = widgets.Checkbox(value=param['is_double_stranded']['default'])
    box = widgets.HBox([widgets.Label(param['is_double_stranded']['glossory'], layout={'width': '400px'}), is_double_stranded])
    display(box)
    options['RuntimeParameters']['is_double_stranded'] =  is_double_stranded



def display_widgets(_options):
    """Display all widgets Jupyter notebook and options

    return widgets which can be used to extract user options
    """

    display(widgets.HTML(value='<b>Backbone</b>'))
    display(widgets.interactive(backbone, param=widgets.fixed(_options['Backbone']),
                                file_path=widgets.Text(value='', description="Backbone File",
                                style={'description_width': 'initial'})))
    display(widgets.HTML(value='<b>Bases</b>'))
    display(widgets.interactive(bases, param=widgets.fixed(_options['Base 1']), 
                                num_bases=widgets.IntSlider(min=1,max=10,step=1,value=1,
                                description="Number of Bases", style={'description_width': 'initial'})))
    display(widgets.HTML(value='<b>Helical Parameters</b>'))
    helical_parameters(_options['HelicalParameters'])
    display(widgets.HTML(value='<b>Runtime Parameters</b>'))
    runtime_parameters(_options['RuntimeParameters'])

    return options


def extract_options(_options):
    """Extract user options from the widgets"""

    user_options = {}
    for k1, val1 in _options.items():
        user_options[k1] = {}
        for k2, val2 in val1.items():
            if isinstance(val2, str):
                user_options[k1][k2] = val2
            elif k1 == 'Backbone':
                user_options[k1][k2] = val2.value

            elif 'Base' in k1:
                user_options[k1][k2] = val2.value

            elif k1 == 'HelicalParameters':
                user_options[k1][k2] = [val2[0].value[0], val2[0].value[1], val2[1].value]

            else:
                if k2 == 'energy_filter':
                    user_options[k1][k2] = [i.value for i in val2]

                else:
                    user_options[k1][k2] = val2.value

    return user_options
