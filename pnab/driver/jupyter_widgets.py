import yaml
import ipywidgets as widgets
from IPython.display import display, Javascript

from pnab.driver.pNAB import pNAB
from pnab.driver import draw
from pnab.driver.options import _options_dict as _options


options = {}


def path(file_path, param):

    import os

    if not os.path.isfile(file_path):
        param['linker']['default'][0] = param['linker']['default'][1] = 1
        param['interconnects']['default'][0] = param['interconnects']['default'][1] = 1
        return

    num_atoms = draw.view_py3dmol(file_path, label=True)

    linker1 = widgets.Dropdown(value=param['linker']['default'][0], options=range(1, num_atoms + 1))
    linker2 = widgets.Dropdown(value=param['linker']['default'][1], options=range(1, num_atoms + 1))
    interconnects1 = widgets.Dropdown(value=param['interconnects']['default'][0], options=range(1, num_atoms + 1))
    interconnects2 = widgets.Dropdown(value=param['interconnects']['default'][1], options=range(1, num_atoms + 1))

    box1 = widgets.HBox([widgets.Label(param['linker']['glossory'], layout={'width': '400px'}), linker1, linker2])
    box2 = widgets.HBox([widgets.Label(param['interconnects']['glossory'],
                layout={'width': '400px'}), interconnects1, interconnects2])

    display(box1)
    display(box2)

    options['Backbone'] = {'file_path': file_path, 'linker': [linker1, linker2],
                           'interconnects': [interconnects1, interconnects2]}



def upload_backbone(f, param):
    if f:
        param['linker']['default'][0] = param['linker']['default'][1] = 1
        param['interconnects']['default'][0] = param['interconnects']['default'][1] = 1

        input_file = list(f.keys())[0]
        with open(input_file, 'wb') as w:
            w.write(list(f.values())[0]['content'])

        param['file_path']['default'] = input_file

    display(widgets.interactive(path, file_path=widgets.Text(value=param['file_path']['default'], description="Backbone File", style={'description_width': 'initial'}),
                                param=widgets.fixed(param), upload=widgets.fixed(True)))


def backbone(param):
    """Backbone widget for use in Jupyter notebook"""

    display(widgets.HTML(value='<b>Backbone</b>'))

    w = widgets.interactive(upload_backbone, param=widgets.fixed(param),
                            f=widgets.FileUpload(accept='', multiple=False, 
                            description="Backbone File", style={'description_width': 'initial'}))
    display(w)


#def base(param, file_path, base_number):
#    """Base widget for use in Jupyter notebook"""
#    import os
#
#    if not os.path.isfile(file_path):
#        return
#
#    num_atoms = draw.view_py3dmol(file_path, label=True)
#
#    linker1 = widgets.Dropdown(options=range(1, num_atoms + 1)) 
#    linker2 = widgets.Dropdown(options=range(1, num_atoms + 1)) 
#    
#    code = widgets.Text()
#    name = widgets.Text()
#    pair_name = widgets.Text()
#
#    box1 = widgets.HBox([widgets.Label(param['linker']['glossory'], layout={'width': '400px'}), linker1, linker2])
#    box2 = widgets.HBox([widgets.Label(param['code']['glossory'], layout={'width': '400px'}), code])
#    box3 = widgets.HBox([widgets.Label(param['name']['glossory'], layout={'width': '400px'}), name])
#    box4 = widgets.HBox([widgets.Label(param['pair_name']['glossory'], layout={'width': '400px'}), pair_name])
#
#    display(box1)
#    display(box2)
#    display(box3)
#    display(box4)
#
#    options['Base %i' %base_number]['file_path'] = file_path
#    options['Base %i' %base_number]['linker'] = [linker1, linker2]
#    options['Base %i' %base_number]['code'] = code
#    options['Base %i' %base_number]['name'] = name
#    options['Base %i' %base_number]['pair_name'] = pair_name
#
#
#def bases(param, num_bases):
#    """Number of bases widget for use in Jupyter notebook
#
#    Dynamically changes the number of base options based on the number of bases
#    """
#
#    for i in range(num_bases):
#        options['Base %i' %(i+1)] = {}
#        path = widgets.interactive(base, param=widgets.fixed(param), file_path=widgets.Text(value='',
#                                   description="Base File", style={'description_width': 'initial'}),
#                                   base_number=widgets.fixed(i+1))
#        display(path)


def bases():
    display(widgets.HTML(value='<b>Bases</b>'))
    display(widgets.HTML(value=('These bases are already defined:' +
                               '<br> Adenine (A), Guanine (G), Cytosine (C), Uracil (U), Thymine (T), ' +
                               'Cyanuric Acid (X), and Triaminopyrimidine (Y)' + 
                               '<br> To add other bases, please edit the generated input yaml file and run it directly.')))



def helical_parameters(param):
    """Helical parameter widget for use in Jupyter notebook"""

    display(widgets.HTML(value='<b>Helical Parameters</b>'))
    display(widgets.HTML(value='All angles are in degrees and all distances in Angstroms.' + 
                               '<br> Random configurations over the ranges of each helical parameter will be generated for sampling.' +
                               '<br> The steps are the number of configurations generated for each helical parameter.'))

    param_dict = {}
    options['HelicalParameters'] = {}
    for k in param:
        param_dict[k] = []
        default = [param[k]['default'][0], param[k]['default'][1], param[k]['default'][2]]
        if k in ['inclination', 'tip', 'twist']:
            param_dict[k].append(widgets.FloatRangeSlider(value=[default[0], default[1]], min=-180, max=180, step=0.01, readout_format='.2f'))
        else:
            param_dict[k].append(widgets.FloatRangeSlider(value=[default[0], default[1]], min=-10, max=10, step=0.01, readout_format='.3f'))
        param_dict[k].append(widgets.BoundedIntText(value=default[2], min=1, max=1000, step=1, description='Steps'))
        box = widgets.HBox([widgets.Label(param[k]['glossory'], layout={'width': '100px'}), param_dict[k][0], param_dict[k][1]])
        display(box)
        options['HelicalParameters'][k] = [param_dict[k][0], param_dict[k][1]]


def runtime_parameters(param):
    """Runtime parameter widget for use in Jupyter notebook"""

    display(widgets.HTML(value='<b>Runtime Parameters</b>'))
    display(widgets.HTML(value='All energies are in kcal/mol and all distances in Angstroms'))

    options['RuntimeParameters'] = {}

    num_steps = widgets.IntText(value=param['num_steps']['default'])
    box = widgets.HBox([widgets.Label(param['num_steps']['glossory'], layout={'width': '400px'}), num_steps])
    display(box)
    options['RuntimeParameters']['num_steps'] = num_steps

    ff_type = widgets.Text(value=param['type']['default'])
    box = widgets.HBox([widgets.Label(param['type']['glossory'], layout={'width': '400px'}), ff_type])
    display(box)
    options['RuntimeParameters']['type'] = ff_type

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

    strand = widgets.Text(value=''.join(param['strand']['default']))
    box = widgets.HBox([widgets.Label(param['strand']['glossory'], layout={'width': '400px'}), strand])
    display(box)
    options['RuntimeParameters']['strand'] = strand


    is_double_stranded = widgets.Checkbox(value=param['is_double_stranded']['default'])
    box = widgets.HBox([widgets.Label(param['is_double_stranded']['glossory'], layout={'width': '400px'}), is_double_stranded])
    display(box)
    options['RuntimeParameters']['is_double_stranded'] =  is_double_stranded

    pair_A_U = widgets.Checkbox(value=param['pair_A_U']['default'])
    box = widgets.HBox([widgets.Label(param['pair_A_U']['glossory'], layout={'width': '400px'}), pair_A_U])
    display(box)
    options['RuntimeParameters']['pair_A_U'] =  pair_A_U

    is_hexad = widgets.Checkbox(value=param['is_hexad']['default'])
    box = widgets.HBox([widgets.Label(param['is_hexad']['glossory'], layout={'width': '400px'}), is_hexad])
    display(box)
    options['RuntimeParameters']['is_hexad'] =  is_hexad

    is_parallel = widgets.Checkbox(value=param['is_parallel']['default'])
    box = widgets.HBox([widgets.Label(param['is_parallel']['glossory'], layout={'width': '400px'}), is_parallel])
    display(box)
    options['RuntimeParameters']['is_parallel'] =  is_parallel


def upload_input(param, f):
    if f:
        input_file = list(f.keys())[0]
        with open(input_file, 'wb') as w:
            w.write(list(f.values())[0]['content'])

        display(widgets.HTML(value=input_file))
        show(param, input_file) 


def show(param, input_file):

    if input_file == "Upload file":
        display(widgets.interactive(upload_input, param=widgets.fixed(param), f=widgets.FileUpload(accept='', multiple=False, description="Input File")))
        return

    options = yaml.load(open(input_file, 'r'), yaml.FullLoader)

    for k1, v1 in options.items():
        for k2, v2 in options[k1].items():
            param[k1][k2]['default'] = options[k1][k2]

    backbone(param['Backbone'])
    bases()
    helical_parameters(param['HelicalParameters'])
    runtime_parameters(param['RuntimeParameters'])


def run(b):
    run_options = extract_options(options)
    with open('options.yaml', 'w') as f:
     f.write(yaml.dump(run_options))

    run = pNAB('options.yaml')
    run.run()
    run.get_results()



def display_widgets(param):
    """Display all widgets Jupyter notebook and options

    return widgets which can be used to extract user options
    """


    disable_js = """
IPython.OutputArea.prototype._should_scroll = function(lines) {
    return false;
}
"""
    display(Javascript(disable_js))
    draw.draw()
    display(widgets.interactive(show, param=widgets.fixed(param),
            input_file=widgets.Dropdown(options=['files/options_rna.yaml', 'files/options_dna.yaml', 'files/options_hexad.yaml', 'Upload file'],
                                        style={'description_width': 'initial'}, description='Input File')))
    button = widgets.Button(description='Run')
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
                if isinstance(val2, list):
                    user_options[k1][k2] = [val2[0].value, val2[1].value]
                else:
                    user_options[k1][k2] = val2

            elif 'Base' in k1:
                if isinstance(val2, list):
                    user_options[k1][k2] = [val2[0].value, val2[1].value]
                else:
                    print(k1, k2, val2)
                    user_options[k1][k2] = val2.value

            elif k1 == 'HelicalParameters':
                user_options[k1][k2] = [val2[0].value[0], val2[0].value[1], val2[1].value]

            else:
                if k2 == 'energy_filter':
                    user_options[k1][k2] = [i.value for i in val2]

                else:
                    user_options[k1][k2] = val2.value

    return user_options


def builder():
    display_widgets(_options)
