"""
This is a file that contains visualization utilities used
with the Jupyter notebook driver. To access these utilities,
internet connection is required. This lowers the number
of dependencies that are required for installation.
The visualization utilities are two: 1) JSME for drawing
molecules. 2) NGLViewer for visualizing results.
Visualization using 3DMol is also included for reference.
"""

import subprocess

from IPython.display import HTML, display

def jsme():
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
        var command = name + ' = ' + `"""` + drawing + `"""`;
        var kernel = IPython.notebook.kernel;
        kernel.execute(command);
    }

    </script>

    <div id="appletContainer"></div>
    <div><button onclick="onSubmit()">done!</button></div>
    '''

    return jsme_html


def process(molecules):
    """
    A function that processes the input molecules drwan by the user.
    A mol file, containing the two-dimensional representation of the molecule,is first written.
    Then, a three-dimensional geometry is created in pdb format using OpenBabel.
    The input geometry is then minimized using OpenBabel. 
    """

    for i in range(len(molecules)):
        name, mol = molecules[i][0].strip(), molecules[i][1]
        f = open(name + '.mol', 'w'); f.write(mol); f.close()
        subprocess.call(['obabel', name + '.mol', '-O', name + '.pdb', '--gen3d', '-d'])
        subprocess.call(['obminimize', name + '.pdb', '>', name + '.pdb'])



def view_3dmol(conformer):
    """
    A function to view molecular geometries using the 3DMol project. Added here for reference only.
    NGLViewer is used here.
    """

    html = """
    <script src="http://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script> </head>    
             <div style="height: 400px; width: 400px; position: relative;" class='viewer_3Dmoljs' data-element='test' data-backgroundcolor='0xffffff' data-style='stick'></div>
    """
    display(HTML("""<textarea id="test" style="display:none;">\n""" + open(conformer).read() + """</textarea>""" + html))


def view_nglviewer(conformer, atom_index='0'):
    """
    A function to view molecular geometries using the NGLViewer project.
    Atom indicies can be added so that the user can pick atoms that connect
    backbones together or connect the bases and the backbones.
    The development version of the Javascript file is used because it has
    functions for showing distances and angles.
    """

    html = """
    <script src="https://unpkg.com/ngl@2.0.0-dev.34/dist/ngl.js"></script>
    <script>
      document.addEventListener("DOMContentLoaded", function () {
        var stage = new NGL.Stage("viewport", { backgroundColor: "white"});
        stage.loadFile("%s", {defaultRepresentation: true}).then(function (component) {
          var lowEnd = 1;
          var highEnd = %s;
          var arr = [];
          while(lowEnd <= highEnd){
                 arr.push(String(lowEnd++));
          }       
    
          component.addRepresentation("label", {labelType:"text", labelText:arr, attachment: 'middle-center', showBorder: true});
        });
      });
    </script>
    <div id="viewport" style="width:500px; height:500px;"></div>
    """ %(conformer, atom_index)

    display(HTML(html, metadata={'isolated': True}))
