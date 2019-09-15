"""!@brief A file containing helper scripts for visualizing molecules

The visualization utilities used with the Jupyter notebook driver is:

1) NGLView for visualizing results.
"""

from __future__ import division, absolute_import, print_function

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
