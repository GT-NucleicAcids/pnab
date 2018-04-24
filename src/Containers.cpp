//
// Created by jbarnett8 on 9/20/17.
//

#include "Containers.h"
#include <openbabel/math/align.h>
#include <openbabel/obconversion.h>

using namespace std;
using namespace PNAB;
using namespace OpenBabel;

void Backbone::deleteVectorAtom() {
    if (!vector_atom_deleted) {
        array<size_t, 3> tmpIDs{backbone.GetAtom(interconnects[0])->GetId(),
                                          backbone.GetAtom(interconnects[1])->GetId(),
                                          backbone.GetAtom(linker[0])->GetId()};
        backbone.DeleteAtom(getVector());
        interconnects = {backbone.GetAtomById(tmpIDs[0])->GetIdx(), backbone.GetAtomById(tmpIDs[1])->GetIdx()};
        linker = {backbone.GetAtomById(tmpIDs[2])->GetIdx(), 0};
        vector_atom_deleted = true;
    }
}

Backbone::Backbone(FileParser &sp) {
    vector_atom_deleted = false;
    Category backbone_params = sp.getCategory("BACKBONE PARAMETERS");
    string file_path = backbone_params.getStringFieldDataFromName("Backbone_File_Path");
    OBConversion conv(file_path);
    OBFormat *inFormat = conv.FormatFromExt(file_path);
    if (inFormat) {
        conv.SetInFormat(inFormat);
    } else {
        cerr << "Backbone: Cannot determine file type from extension. Defaulting to PDB." << endl;
        conv.SetInFormat("PDB");
    }
    if (!conv.Read(&backbone)) {
        cerr << "Backbone: There was an error reading file for backbone: " << file_path << endl;
        exit(1);
    }

    auto vec = backbone_params.getSizeVecFieldDataFromName("Interconnects");
    if (vec.size() != 2) {
        cerr << "Backbone: Incorrect number of elements specified for Field \"Interconnects\" in Category \""
             << backbone_params.getName() << "\". There must be exactly 2 indices indicated." << endl;
        exit(1);
    }
    copy(vec.begin(), vec.end(), interconnects.begin());

    vec = backbone_params.getSizeVecFieldDataFromName("Base_Connect");
    if (vec.size() != 2) {
        cerr << "Backbone: Incorrect number of elements specified for Field \"Base_Connect\" in Category \""
             << backbone_params.getName() << "\". There must be exactly 2 indices indicated." << endl;
        exit(1);
    }
    copy(vec.begin(), vec.end(), linker.begin());
    validate();
}

void Backbone::validate() {
    if (interconnects.size() == 2 & linker.size() == 2) {
        for (auto i : interconnects)
            if (i > backbone.NumAtoms()) {
                cerr << "Backbone: Interconnects atom indices are out of bounds." << endl;
                exit(1);
            }

        for (auto i : linker)
            if (i > backbone.NumAtoms()) {
                cerr << "Backbone: Linker atom indices are out of bounds." << endl;
                exit(1);
            }

        OBAtom *link = backbone.GetAtom(static_cast<int>(linker[0])),
                *vec = backbone.GetAtom(static_cast<int>(linker[1]));
        if (!link->IsConnected(vec)) {
            cerr << "Backbone: Atoms in defined in the Linker must be neighbors." << endl;
            exit(1);
        }

    } else {
        cerr << "Backbone: Interconnects or Linker is of incorrect dimension." << endl;
        exit(1);
    }
}

void Base::validate() {
    if (linker.size() == 2) {
        for (auto i : linker)
            if (i > base.NumAtoms()) {
                cerr << "Base: Linker atom indices are out of bounds." << endl;
                exit(1);
            }

        OBAtom *link = base.GetAtom(static_cast<int>(linker[0])),
               *vec = base.GetAtom(static_cast<int>(linker[1]));
        if (!link->IsConnected(vec)) {
            cerr << "Base: Atoms in defined in the Linker must be neighbors." << endl;
            exit(1);
        }

    } else {
        cerr << "Base: Linker is of incorrect dimension." << endl;
        exit(1);
    }
}

Bases::Bases(FileParser &fp) {
    Category base_params = fp.getCategory("BASE PARAMETERS");
    vector<string> file_paths = base_params.getStringVecFieldDataFromName("Base_File_Path");
    if (file_paths.size() == 0) {
        cerr << "Bases: No file paths were specified." << endl;
        exit(1);
    }
    size_t num_bases = file_paths.size(), i = 0;
    vector<OBMol> mols(num_bases);
    OBConversion conv;
    for (auto file : file_paths) {
        OBFormat *inFormat = conv.FormatFromExt(file);
        if (inFormat) {
            conv.SetInFormat(inFormat);
        } else {
            cerr << "Base: Cannot determine file type from extension. Defaulting to PDB." << endl;
            conv.SetInFormat("PDB");
        }
        if (!conv.ReadFile(&mols[i], file)) {
            cerr << "Base: There was an error reading file for backbone: " << file << endl;
            exit(1);
        }
        i++;
    }

    auto vec = base_params.getSizeVecFieldDataFromName("Backbone_Connect");
    if (vec.size() != 2 * num_bases) {
        cerr << "Incorrect number of elements specified for Field \"Backbone_Connect\" in Category \""
             << base_params.getName() << "\". There must be exactly 2 * (number of bases) indices indicated." << endl;
        exit(1);
    }

    auto name_vec = base_params.getStringVecFieldDataFromName("Name");
    auto code_vec = base_params.getStringVecFieldDataFromName("Code");
    if (name_vec.size() != num_bases) {
        cerr << "Incorrect number of names specified in Category 'BASE PARAMETERS'." << endl;
        exit(1);
    }
    if (code_vec.size() != num_bases) {
        cerr << "Incorrect number of codes specified in Category 'BASE PARAMETERS'." << endl;
        exit(1);
    }

    for (size_t j = 0; j < num_bases; ++j) {
        array<size_t, 2> two{vec[2*j], vec[2*j + 1]};
        transform(code_vec[j].begin(), code_vec[j].end(), code_vec[j].begin(), ::toupper);
        bases.push_back(Base(name_vec[j], code_vec[j], mols[j], two));
    }

}

BaseUnit::BaseUnit(Base base, Backbone backbone) {

    // We go ahead and center everything so we can rotate later
    backbone.center();

    // Setting up an array of necessary \code{OBAtom}s for alignment
    array< OBAtom*, 4 > atoms{base.getLinker() , base.getVector()  ,
                                   backbone.getLinker(), backbone.getVector() };

    vector<vector3> ref {atoms[1]->GetVector(), atoms[0]->GetVector()},
            target {atoms[2]->GetVector(), atoms[3]->GetVector()},
            result;

    // Aligning backbone linker to base linker
    OBAlign align(ref, target);
    matrix3x3 matrix;
    if (align.Align()) {
        matrix = align.GetRotMatrix();
    } else {
        cerr << "Failed to align backbone to base. Check base and backbone." << endl;
        exit(1);
    }
    double rot[9];
    matrix.GetArray(rot);


    // garbage code for debugging
    // TODO delete this code
//    {
//        OBConversion conv;
//        std::filebuf fb;
//        fb.open("center_out.pdb", std::ios::out);
//        std::ostream fileStream(&fb);
//        conv.SetOutFormat("PDB");
//        conv.SetOutStream(&fileStream);
//        OBMol mol;
//        // mol += base.getMolecule();
//        mol += backbone.getMolecule();
//        conv.Write(&mol);
//    }

    backbone.rotate(rot);
    backbone.translate(atoms[0]->GetVector() - atoms[3]->GetVector());

    // garbage code for debugging
    // TODO delete this code
//    {
//        OBConversion conv;
//        std::filebuf fb;
//        fb.open("rot_trans_out.pdb", std::ios::out);
//        std::ostream fileStream(&fb);
//        conv.SetOutFormat("PDB");
//        conv.SetOutStream(&fileStream);
//        OBMol mol;
//        mol += base.getMolecule();
//        mol += backbone.getMolecule();
//        conv.Write(&mol);
//    }

    // Deleting old hydrogens that formed the vector
    backbone.deleteVectorAtom();
    base.deleteVectorAtom();

    // garbage code for debugging
    // TODO delete this code
//    {
//        OBConversion conv;
//        std::filebuf fb;
//        fb.open("rot_trans_delete_atoms_out.pdb", std::ios::out);
//        std::ostream fileStream(&fb);
//        conv.SetOutFormat("PDB");
//        conv.SetOutStream(&fileStream);
//        OBMol mol;
//        mol += base.getMolecule();
//        mol += backbone.getMolecule();
//        conv.Write(&mol);
//    }

    OBMol mol;
    mol += base.getMolecule();

    base_bond_indices = {1, mol.NumBonds()};

    unsigned num_atoms = mol.NumAtoms();
    mol += backbone.getMolecule();
    mol.AddBond(base.getLinker()->GetIdx(), backbone.getLinker()->GetIdx() + num_atoms, 1);

    // garbage code for debugging
    // TODO delete this code
    {
        OBConversion conv;
        std::filebuf fb;
        fb.open("base_unit_example_out.pdb", std::ios::out);
        std::ostream fileStream(&fb);
        conv.SetOutFormat("PDB");
        conv.SetOutStream(&fileStream);
        OBMol mol;
        mol += base.getMolecule();
        mol += backbone.getMolecule();
        mol.AddBond(base.getLinker()->GetIdx(), backbone.getLinker()->GetIdx() + num_atoms, 1);

        OBAtom *a = mol.GetAtom(10), *b = mol.GetAtom(5), *c = mol.GetAtom(22), *d = mol.GetAtom(18);
        mol.SetTorsion(a, b, c, d, 82.2 * DEG_TO_RAD);

        conv.Write(&mol);
    }

    unit = mol;
    base_connect_index = base.getLinker()->GetIdx();
    backbone_interconnects = {backbone.getHead()->GetIdx() + num_atoms, backbone.getTail()->GetIdx() + num_atoms};
    base_index_range = {1, num_atoms};
    backbone_index_range = {num_atoms + 1, mol.NumAtoms()};
}

RuntimeParameters::RuntimeParameters(FileParser &sp) {
    Category rp = sp.getCategory("RUNTIME PARAMETERS");

    // Geometric Parameters
    base_to_backbone_bond_length = rp.getDoubleFieldDataFromName("Base_to_Backbone_Bond_Length", numeric_limits<double>::quiet_NaN());

    // Energetic parameters
    vector<string> single_doubles_energy{"Max_Total_Energy", "Max_Angle_Energy", "Max_Bond_Energy",
                                         "Max_VDW_Energy", "Max_Torsion_Energy"};
    energy_filter.clear();
    for (auto d_val : single_doubles_energy)
        energy_filter.push_back(rp.getDoubleFieldDataFromName(d_val));
    max_distance = rp.getDoubleFieldDataFromName("Max_Backbone_Interlink_Distance");

    // ForceField Parameters
    type = rp.getStringFieldDataFromName("Force_Field_Type","GAFF");
    parameter_file = rp.getStringFieldDataFromName("Force_Field_Parameter_File");
    num_steps = rp.getSizeFieldDataFromName("Search_Size");
    dihedral_discretization = rp.getSizeFieldDataFromName("Dihedral_Step_Size",numeric_limits<size_t>::quiet_NaN());
    angleStepSize = rp.getSizeFieldDataFromName("Search_Step_Size",numeric_limits<size_t>::quiet_NaN());
    chain_length = rp.getSizeFieldDataFromName("Chain_Length",3ul);
    algorithm = rp.getStringFieldDataFromName("Algorithm");
    validate();
}

void RuntimeParameters::validate() {

    if (!isnan(base_to_backbone_bond_length) & base_to_backbone_bond_length < 0) {
        cerr << "Base_to_Backbone_Bond_Length must be greater than 0." << endl;
        exit(1);
    }
    if (max_distance < 0) {
        cerr << "Max_Distance must be greater than 0." << endl;
        exit(1);
    }
}

HelicalParameters::HelicalParameters(FileParser &fp) {
    Category hp = fp.getCategory("HELICAL PARAMETERS");

    vector<RangeField*> ranges{&tilt, &roll, &twist, &shift, &slide, &rise, &buckle, &propeller,
                               &opening, &shear, &stretch, &stagger, &inclination, &tip,
                               &x_displacement, &y_displacement};
    vector<string> names{"Tilt", "Roll", "Twist", "Shift", "Slide", "Rise", "Buckle", "Propeller", "Opening",
                         "Shear", "Stretch", "Stagger", "Inclination", "Tip", "X_Displacement", "Y_Displacement"};

    unsigned i = 0;
    for (auto &r : ranges) {
        r->d = hp.getDoubleVecFieldDataFromName(names[i++]);
        checkRangeField(*r);
    }
}
