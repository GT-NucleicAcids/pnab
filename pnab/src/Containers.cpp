#include "Containers.h"
#include <openbabel/math/align.h>

using namespace std;
using namespace PNAB;
using namespace OpenBabel;

void Backbone::deleteVectorAtom() {
    // Delete the extra atom that defines the vector connecting the backbone to the base
    // If it is already deleted, do nothing
    if (!vector_atom_deleted) {
        // Save the IDs of atoms connecting the two backbones and the ID of the atom connecting to the base
        array<size_t, 3> tmpIDs{backbone.GetAtom(interconnects[0])->GetId(),
                                          backbone.GetAtom(interconnects[1])->GetId(),
                                          backbone.GetAtom(linker[0])->GetId()};
        // Delete the atom
        backbone.DeleteAtom(getVector());
        // Get the new indices after the deletion of the atom
        interconnects = {backbone.GetAtomById(tmpIDs[0])->GetIdx(), backbone.GetAtomById(tmpIDs[1])->GetIdx()};
        linker = {backbone.GetAtomById(tmpIDs[2])->GetIdx(), 0};
        vector_atom_deleted = true;
    }
}

Backbone::Backbone(std::string file_path, std::array<unsigned, 2> interconnects, std::array<unsigned,2> linker, std::vector<std::vector<unsigned>> fixed_bonds) {
    // Constructor for the backbone

    vector_atom_deleted = false;
    this->interconnects = interconnects;
    this->linker = linker;
    this->fixed_bonds = fixed_bonds;

    // Read backbone file
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

    // Basic validation of the indices
    validate();
}

void Backbone::validate() {

    // Make sure only two atoms are specified for connecting the backbones
    auto vec = interconnects;
    if (vec.size() != 2) {
        cerr << "Backbone: Incorrect number of elements specified for the atoms connecting the two backbones."
             << " There must be exactly 2 indices indicated." << endl;
        exit(1);
    }
    copy(vec.begin(), vec.end(), interconnects.begin());

    // Make sure only two atoms are specified for connecting backbone to base
    vec = linker;
    if (vec.size() != 2) {
        cerr << "Backbone: Incorrect number of elements specified for the atoms forming the vector connecting the backbone to the nucleobase."
             << " There must be exactly 2 indices indicated." << endl;
        exit(1);
    }
    copy(vec.begin(), vec.end(), linker.begin());

    for (auto i : interconnects)
        // Make sure the index is within the number of atoms in the backbone
        if (i > backbone.NumAtoms()) {
            cerr << "Backbone: Interconnects atom indices are out of bounds." << endl;
            exit(1);
        }

    for (auto i : linker)
        // Make sure the index is within the number of atoms in the backbone
        if (i > backbone.NumAtoms()) {
            cerr << "Backbone: Linker atom indices are out of bounds." << endl;
            exit(1);
        }

    OBAtom *link = backbone.GetAtom(static_cast<int>(linker[0])),
            *link2 = backbone.GetAtom(static_cast<int>(linker[1]));

    // Make sure the two provided atoms are neighbors
    if (!link->IsConnected(link2)) {
        cerr << "Backbone: Atoms in defined in the Linker must be neighbors." << endl;
        exit(1);
    }

}


Base::Base(std::string nameT, std::string codeT, std::string file_path, std::array<std::size_t, 2> linkerT, std::string pairT) {
    // Constructor for the base

    name = nameT;
    code = codeT;
    linker = linkerT;

    if (file_path.empty()) {
        cerr << "Bases: No file path was specified." << endl;
        exit(1);
    }

    // Read file specifying the geometry of the base
    OBConversion conv;
    OBFormat *in_format = conv.FormatFromExt(file_path);
    if (in_format)
        conv.SetInFormat(in_format);
    else {
        cerr << "Base: Cannot determine file type from extension. Defaulting to PDB." << endl;
        conv.SetInFormat("PDB");
    }
    OBMol molT;
    if (!conv.ReadFile(&molT, file_path)) {
        cerr << "Base: There was an error reading file for base: " << file_path << endl;
        exit(1);
    }

    base = OBMol(molT);
    vector_atom_deleted = false;
    pair_name = pairT;

    // validate
    validate();
}


void Base::validate() {
    // Basic validation

    // Make sure we have the correct size for the atom indices
    if (linker.size() == 2) {
        for (auto i : linker)
            // Make sure the indices are within the number of atoms in the base
            if (i > base.NumAtoms()) {
                cerr << "Base: Linker atom indices are out of bounds." << endl;
                exit(1);
            }

        OBAtom *link = base.GetAtom(static_cast<int>(linker[0])),
               *vec = base.GetAtom(static_cast<int>(linker[1]));

        // Make sure the atoms are neighbors
        if (!link->IsConnected(vec)) {
            cerr << "Base: Atoms in defined in the Linker must be neighbors." << endl;
            exit(1);
        }

    } else {
        cerr << "Base: Linker is of incorrect dimension." << endl;
        exit(1);
    }
}

Bases::Bases(std::vector<Base> input_bases) {

    // Loop through the bases and create a vector of them 
    for (unsigned i = 0; i < input_bases.size(); ++i) {
        string file_path = input_bases[i].file_path;

        auto vec = input_bases[i].linker;
        if (vec.size() != 2) {
            cerr << "Incorrect number of elements specified for Field \"Backbone_Connect\" "
                 << "\". There must be exactly 2 indices indicated." << endl;
            exit(1);
        }
        array<size_t, 2> linkers = {vec[0], vec[1]};

        auto name = input_bases[i].name;
        transform(name.begin(), name.end(), name.begin(), ::tolower);
        auto code = input_bases[i].code;
        transform(code.begin(), code.end(), code.begin(), ::toupper);
        auto pair_name = input_bases[i].pair_name;
        transform(pair_name.begin(), pair_name.end(), pair_name.begin(), ::tolower);

        // Call Base again to make sure we have the openbabel molecules
        bases.push_back(Base(name,code,file_path,linkers, pair_name));
    }

    // Checking to see if the bases have pairs
    all_bases_pair = true;
    for (auto v : bases) {
        auto pair_name = v.getBasePairName();
        transform(pair_name.begin(), pair_name.end(), pair_name.begin(), ::tolower);
        bool pair_found = false;
        for (auto w : bases) {
            if (w.getName().find(pair_name) != string::npos) {
                auto base_name = v.getName();
                transform(base_name.begin(), base_name.end(), base_name.begin(), ::tolower);
                name_base_map.insert(pair<string,Base>(v.getName(),w));
                pair_found = true;
                break;
            }
        }
        if (!pair_found) {
            all_bases_pair = false;
            break;
        }
    }
}

std::vector<Base> Bases::getBasesFromStrand(std::vector<std::string> strand) {
    // Return the bases given the names of the bases in the strand
    vector<Base> bases;
    for (auto &v : strand) {
        transform(v.begin(), v.end(), v.begin(), ::tolower);
        bases.push_back(getBaseFromName(v));
    }
    return bases;
}

std::vector<Base> Bases::getComplimentBasesFromStrand(std::vector<std::string> strand) {
    // Return the complimentary bases given the names of the bases in the strand
    if (all_bases_pair) {
        vector<Base> compliment_bases;
        for (auto &v : strand) {
            transform(v.begin(), v.end(), v.begin(), ::tolower);
            compliment_bases.push_back(name_base_map.find(v)->second);
        }
        return compliment_bases;
    }

    // If not all bases have complimentary bases, return empty vector
    return vector<Base>();
}

BaseUnit::BaseUnit(Base base, Backbone backbone) {

    // We go ahead and center everything so we can rotate later
    backbone.center();

    // Setting up an array of necessary OBAtom for alignment
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

    backbone.rotate(rot);

    // Necessary to invert to get correct DNA and RNA structures
    // if alignment requires improper rotation
    if (matrix.determinant() < 0) {
        double invert[9] = {1, 0, 0, 0, 1, 0, 0, 0, -1};
        backbone.rotate(invert);
    }

    backbone.translate(atoms[0]->GetVector() - atoms[3]->GetVector());

    // Get indices of fixed atoms
    unsigned num_atoms = base.getMolecule().NumAtoms() - 1;
    for (int i=0; i < backbone.fixed_bonds.size(); i++) {
        fixed_bonds.push_back({});
        for (auto j: backbone.fixed_bonds[i]) {
            if (j == backbone.linker[1]) {
                // This atom will be deleted, so put the index of the atom from the base
                if (base.linker[0] < base.linker[1])
                    fixed_bonds[i].push_back(base.linker[0]);
                else
                    fixed_bonds[i].push_back(base.linker[0]-1);
            }
            else if (j < backbone.linker[1])
                // The index of the atom will not change because of deletion
                fixed_bonds[i].push_back(j + num_atoms);
            else
                // Subtract one because of the deletion of the atom
                fixed_bonds[i].push_back(j-1 + num_atoms);
        }
    }

    // Deleting old atoms that formed the vector
    backbone.deleteVectorAtom();
    base.deleteVectorAtom();

    // Create the molecule
    OBMol mol;
    mol += base.getMolecule();

    num_atoms = mol.NumAtoms();
    mol += backbone.getMolecule();

    // Form a bond between base and backbone
    mol.AddBond(base.getLinker()->GetIdx(), backbone.getLinker()->GetIdx() + num_atoms, 1);

    // Set length of new bond to equilibrium length (the sum of the van der Waals radii)
    OBBond* new_bond = mol.GetBond(base.getLinker()->GetIdx(), backbone.getLinker()->GetIdx() + num_atoms);
    new_bond->SetLength(mol.GetAtom(base.getLinker()->GetIdx()), new_bond->GetEquibLength());

    // garbage code for debugging
    // TODO delete this code
    //{
    //    OBConversion conv;
    //    std::filebuf fb;
    //    fb.open("base_unit_example_out.pdb", std::ios::out);
    //    std::ostream fileStream(&fb);
    //    conv.SetOutFormat("PDB");
    //    conv.SetOutStream(&fileStream);
    //    OBMol mol;
    //    mol += base.getMolecule();
    //    mol += backbone.getMolecule();
    //    mol.AddBond(base.getLinker()->GetIdx(), backbone.getLinker()->GetIdx() + num_atoms, 1);
    //    OBBond* new_bond = mol.GetBond(base.getLinker()->GetIdx(), backbone.getLinker()->GetIdx() + num_atoms);
    //    new_bond->SetLength(mol.GetAtom(base.getLinker()->GetIdx()), new_bond->GetEquibLength());
    //    conv.Write(&mol);
    //}

    // Set variables
    unit = mol;
    base_connect_index = num_atoms + backbone.getLinker()->GetIdx();
    backbone_interconnects = {backbone.getHead()->GetIdx() + num_atoms, backbone.getTail()->GetIdx() + num_atoms};
    base_index_range = {1, num_atoms};
    backbone_index_range = {num_atoms + 1, mol.NumAtoms()};
}
