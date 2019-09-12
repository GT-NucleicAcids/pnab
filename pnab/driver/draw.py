"""!@brief A file containing helper scripts for drawing and visualizing molecules

The visualization utilities used with the Jupyter notebook driver are:

1) JSME for drawing molecules\n
2) NGLView for visualizing results.
"""

from __future__ import division, absolute_import, print_function

def draw():
    """!@brief Display the JSME molecule editor (https://peter-ertl.com/jsme/).
    
    HTML script for displaying the JSME molecule editor in the Jupyter notebook.
    After the user draws a molecule and presses 'done!', a prompt asks for the
    name of the molecule. The drawn molecule is processed by OpenBabel to make
    a 3D geometry of the moleucle. The geometry of the molecule is optimized by
    performing 500 steepest descent steps, followed by weighted rotor search,
    followed by 500 conjugate gradient steps. The MMFF94 force field is used
    for the optimization. After the optimization, a PDB file with the user-provided
    name and containing the optimized geometry is written.

    @attention The JSME editor may not be displayed in some browsers.

    @sa jupyter_widgets.backbone
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

    import os
    from IPython.display import display

    import openbabel
    import nglview

    # Check if the file exists
    if not os.path.isfile(molecule):
        return 0

    view = nglview.NGLWidget()

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
        view.add_component(mol, ext='pdb', defaultRepresentation=False) 
        view.add_representation('licorice')
        view.center()
        view.add_representation('label', labelType='serial', backgroundColor='black', showBackground=True)
        display(view)

        # Return the number of atoms to be used for the backbone widgets
        return num_atoms

    # Show accepted candidates
    else:
        # The C++ code always generates PDB files
        view.add_component(molecule, ext='pdb', defaultRepresentation=False)
        view.add_representation('licorice')
        view.center()
        display(view)
