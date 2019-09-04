"""visualization support
Provides visualization utilities.

This is a file that contains visualization utilities used with the Jupyter notebook driver.
The visualization utilities are two: 1) JSME for drawing
molecules. 2) Py3DMol for visualizing results.
"""

from __future__ import division, absolute_import, print_function

def draw():
    """JSME molecule editor.
    
    HTML script for JSME. This has to be displayed from inside the 
    Jupyter notebook in order for the program to get the molecules
    drawn by the user. Jupyter implements security measures that prevents
    the communication of user input to the python driver unless this script
    is called from within the Jupyter notebook. After the user draws a molecule
    and presses 'done!', a prompt asks for the name of the molecule, which becomes
    the variable that holds the mol representation of the geometry. 
    """

    from IPython.display import HTML, display


    jsme_html = '''
    <script type="text/javascript" src="https://peter-ertl.com/jsme/JSME_2017-02-26/jsme/jsme.nocache.js"></script>

    <script>

    function jsmeOnLoad() {
        //arguments: HTML id, width, height (must be string not number!)

        document.JME = new JSApplet.JSME("appletContainer", "380px", "340px", {
                         //optional parameters
                         "options" : "query,noreaction,paste"
        });

    }

    function onSubmit() {
        var drawing = document.JME.molFile();
        var name = prompt("Please enter the name of the molecule. A PDB file with the given name will be generated.", "mol");
        var command = `try: molecules["""` + name + `"""] = """` + drawing + `"""` + `\nexcept: molecules = {"""` + name + `""": """` + drawing + `"""}`; 
        var kernel = IPython.notebook.kernel;
        kernel.execute(command);

        var command = `
import openbabel
mol = openbabel.OBMol()
conv = openbabel.OBConversion()
conv.SetInAndOutFormats("mol", "pdb")
conv.ReadString(mol, """` + drawing + `""")
builder = openbabel.OBBuilder()
builder.Build(mol)
mol.AddHydrogens()
ff = openbabel.OBForceField.FindForceField("mmff94")
ff.Setup(mol)
ff.SteepestDescent(500, 1e-4)
ff.WeightedRotorSearch(250, 10)
ff.ConjugateGradients(500, 1e-6)
ff.GetCoordinates(mol)
with open("""` + name + `""" + '.pdb', 'w') as f:
    f.write(conv.WriteString(mol))
        `;
        var kernel = IPython.notebook.kernel;
        kernel.execute(command);

    }

    </script>

    <div id="appletContainer"></div>

    
    <div><button onclick="onSubmit()">done!</button></div>
    '''


    display(HTML(jsme_html))


def view_nglviwe(conformer, label=False):
    """NGLView viewer.

    A function to view molecular geometries using the NGLView project.
    """

    try:
        import nglview
    except ImportError:
        return

    from IPython.display import display
    import openbabel
    import os
    if not os.path.isfile(conformer):
        return 0

    mol = openbabel.OBMol()
    conv = openbabel.OBConversion()
    fmt = conv.FormatFromExt(conformer)
    conv.SetInAndOutFormats(fmt.GetID(), 'pdb')
    conv.ReadFile(mol, conformer)
    num_atoms = mol.NumAtoms()
    mol = conv.WriteString(mol)

    view = nglview.NGLWidget()
    view.add_component(mol, ext='pdb', defaultRepresentation=False) 
    view.add_representation('licorice')
    if label:
        view.add_representation('label', labelType='serial', backgroundColor='black', showBackground=True)
    view.center()
    display(view)

    return num_atoms
