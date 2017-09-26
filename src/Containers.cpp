//
// Created by jbarnett8 on 9/20/17.
//

#include "Containers.h"
#include <openbabel/math/align.h>

using namespace std;
using namespace PNAB;
using namespace OpenBabel;

Backbone::Backbone(OBMol &molP, std::array<unsigned, 2> interconnectsT, array<unsigned, 2> linkerT) {
    backbone = OBMol(molP);
    interconnects = interconnectsT;
    linker = linkerT;
    vector_atom_deleted = false;
}

OBAtom* Backbone::getHead() {
    return backbone.GetAtom(interconnects[0]);

}
OBAtom* Backbone::getTail() {
    return backbone.GetAtom(interconnects[1]);
}

OBAtom* Backbone::getLinker() {
    return backbone.GetAtom(linker[0]);
}

OBAtom* Backbone::getVector() {
    if (!vector_atom_deleted)
        return backbone.GetAtom(linker[1]);
    else
        return nullptr;
}

void Backbone::center() { backbone.Center(); }

void Backbone::rotate(double* rot) { backbone.Rotate(rot); }

void Backbone::translate(vector3 vec) { backbone.Translate(vec); }

OBMol Backbone::getMolecule() { return backbone; }

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

Base::Base(string nameT, string codeT, OBMol &molT, array<unsigned, 2> linkerT) {
    name = nameT;
    code = codeT;
    linker = linkerT;
    base = OBMol(molT);
    vector_atom_deleted = false;
}

OBAtom* Base::getLinker() {
    return base.GetAtom(linker[0]);
}

OBAtom* Base::getVector() {
    if (!vector_atom_deleted)
        return base.GetAtom(linker[1]);
    else
        return nullptr;
}

OBMol Base::getMolecule() {
    return base;
}

void Base::deleteVectorAtom() {
    if (!vector_atom_deleted) {
        size_t id = getLinker()->GetId();
        base.DeleteAtom(getVector());
        linker = {base.GetAtomById(id)->GetIdx(), 0};
        vector_atom_deleted = true;
    }
}

string Base::getCode() {
    return code;
}

string Base::getName() {
    return name;
}

BaseUnit::BaseUnit(Base b, Backbone bb) {
    // Setting up an array of necessary \code{OBAtom}s for alignment
    std::array< OBAtom*, 4 > atoms{b.getLinker() , b.getVector()  ,
                                   bb.getLinker(), bb.getVector() };
    vector<vector3> ref    {atoms[1]->GetVector(), atoms[0]->GetVector()},
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

    // Rotating and translating backbone to line up with base
    bb.center();
    bb.rotate(rot);
    bb.translate(atoms[0]->GetVector() - atoms[1]->GetVector());

    // Deleting old hydrogens that formed the vector
    bb.deleteVectorAtom();
    b.deleteVectorAtom();

    OBMol mol(b.getMolecule());
    unsigned num_atoms = mol.NumAtoms();
    mol += bb.getMolecule();
    mol.AddBond(b.getLinker()->GetIdx(), bb.getLinker()->GetIdx(), 1);

    unit = OpenBabel::OBMol(mol);
    baseAtomBegin     = mol.GetFirstAtom();
    baseAtomEnd       = mol.GetAtom(b.getVector()->GetIdx());
    backboneAtomBegin = mol.GetAtom(mol.NumAtoms() - num_atoms);
    backboneAtomEnd   = mol.GetAtom(mol.NumAtoms());
}
