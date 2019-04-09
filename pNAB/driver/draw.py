"""
This is a file that contains visualization utilities used
with the Jupyter notebook driver.
The visualization utilities are two: 1) JSME for drawing
molecules. 2) Py3DMol for visualizing results.
"""

from __future__ import division, absolute_import, print_function

def draw():
    """
    HTML script for JSME. This has to be displayed from inside the 
    Jupyter notebook in order for the program to get the molecules
    drawn by the user. Jupyter implements security measures that prevents
    the communication of user input to the python driver unless this script
    is called from within the Jupyter notebook. After the user draws a molecule
    and presses 'done!', a prompt asks for the name of the molecule, which becomes
    the variable that holds the mol representation of the geometry. 
    """

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
        var name = prompt("Please enter the name of the molecule", "mol");
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
ff.SteepestDescent(500)
ff.GetCoordinates(mol)
mol.DeleteHydrogens()
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


    # JSME with highlight
    jsme_html2 = '''

    <script type="text/javascript" src="https://peter-ertl.com/jsme/JSME_2017-02-26/jsme/jsme.nocache.js"></script>

    <script>


    function atomHighLight(molIndex, atomIndex) {
       //index must start at 1
        document.getElementById("atomHighLightTextAreaOut").value = atomIndex;
    }

    function jsmeOnLoad() {
        //arguments: HTML id, width, height (must be string not number!)

        document.JME = new JSApplet.JSME("appletContainer", "380px", "340px", {
                         //optional parameters
                         "options" : "query,noreaction,paste"
        });

        document.JME.setNotifyAtomHighLightChangeJSfunction("atomHighLight")

    }

    function onSubmit() {
        var drawing = document.JME.molFile();
        var name = prompt("Please enter the name of the molecule", "mol");
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
ff.SteepestDescent(1000)
ff.GetCoordinates(mol)
mol.DeleteHydrogens()
with open("""` + name + `""" + '.pdb', 'w') as f:
    f.write(conv.WriteString(mol))
        `;
        var kernel = IPython.notebook.kernel;
        kernel.execute(command);

    }

    </script>

    <body>
    <table><tbody><tr>
    <td id="appletContainer"></td>
    </tr></tbody></table>
    
    <table><tbody><tr>
    <td>Atom highlighted: <textarea id="atomHighLightTextAreaOut" rows="1" cols="3"></textarea></td>
    </tr></tbody></table>
    </body>
    
    <div><button onclick="onSubmit()">done!</button></div>

    '''

    from IPython.display import HTML, display

    display(HTML(jsme_html))


def view_py3dmol(conformer, label=False):
    """
    A function to view molecular geometries using the Py3DMol project.
    """

    try:
        import py3Dmol
    except ImportError:
        return

    import openbabel

    mol = openbabel.OBMol()
    conv = openbabel.OBConversion()
    fmt = conv.FormatFromExt(conformer)
    conv.SetInAndOutFormats(fmt.GetID(), 'pdb')
    conv.ReadFile(mol, conformer)
    num_atoms = mol.NumAtoms()
    mol = conv.WriteString(mol)

    view = py3Dmol.view()
    view.addModel(mol, 'pdb')
    view.setStyle({'stick':{}})
    view.zoomTo()

    if label:
        for i in range(num_atoms):
            view.addLabel(str(i + 1), {}, {'serial': i + 1})

    view.show()


#def view_nglview(conformer, num_atoms=0):
#    """
#    A function to view molecular geometries using the NGLView project.
#    """
#
#    import nglview
#
#    view = nglview.show_file(conformer)
#
#    if num_atoms:
#        label = [str(i) for i in range(1, num_atoms + 1)]
#        view.add_label(selection='all', labelType="text", labelText=label, attachment='middle-center', showBorder=True)
#
#    view.display()
#
#
#def view_3dmol(conformer):
#    """
#    A function to view molecular geometries using the 3DMol project. Added here for reference only.
#    NGLViewer is used here.
#    """
#
#    html = """
#    <script src="http://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script> </head>    
#             <div style="height: 400px; width: 400px; position: relative;" class='viewer_3Dmoljs' data-element='test' data-backgroundcolor='0xffffff' data-style='stick'></div>
#    """
#    display(HTML("""<textarea id="test" style="display:none;">\n""" + open(conformer).read() + """</textarea>""" + html))
#
#
#def view_nglviewer(conformer, atom_index='0'):
#    """
#    A function to view molecular geometries using the NGLViewer project.
#    Atom indicies can be added so that the user can pick atoms that connect
#    backbones together or connect the bases and the backbones.
#    The development version of the Javascript file is used because it has
#    functions for showing distances and angles.
#    """
#
#    html = """
#    <script src="https://unpkg.com/ngl@2.0.0-dev.34/dist/ngl.js"></script>
#    <script>
#      document.addEventListener("DOMContentLoaded", function () {
#        var stage = new NGL.Stage("viewport", { backgroundColor: "white"});
#        stage.loadFile("%s", {defaultRepresentation: true}).then(function (component) {
#          var lowEnd = 1;
#          var highEnd = %s;
#          var arr = [];
#          while(lowEnd <= highEnd){
#                 arr.push(String(lowEnd++));
#          }       
#    
#          component.addRepresentation("label", {labelType:"text", labelText:arr, attachment: 'middle-center', showBorder: true});
#        });
#      });
#    </script>
#    <div id="viewport" style="width:500px; height:500px;"></div>
#    """ %(conformer, atom_index)
#
#    display(HTML(html, metadata={'isolated': True}))
