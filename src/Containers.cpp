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

Backbone::Backbone(std::string file_path, std::array<unsigned, 2> interconnects, std::array<unsigned,2> linker) {
    vector_atom_deleted = false;
    this->interconnects = interconnects;
    this->linker = linker;
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

    auto vec = interconnects;
    if (vec.size() != 2) {
        cerr << "Backbone: Incorrect number of elements specified for Field \"Interconnects\" \""
             << "\". There must be exactly 2 indices indicated." << endl;
        exit(1);
    }
    copy(vec.begin(), vec.end(), interconnects.begin());

    vec = linker;
    if (vec.size() != 2) {
        cerr << "Backbone: Incorrect number of elements specified for Field \"Base_Connect\""
             << "\". There must be exactly 2 indices indicated." << endl;
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


Base::Base(string nameT, string codeT, string file_path, array<std::size_t, 2> linkerT, string pairT) {
    name = nameT;
    code = codeT;
    linker = linkerT;

    if (file_path.empty()) {
        cerr << "Bases: No file path was specified." << endl;
        exit(1);
    }
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
    validate();
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

Bases::Bases(vector<Base> py_bases) {

    for (unsigned i = 0; i < py_bases.size(); ++i) {
        string file_path = py_bases[i].file_path;

        auto vec = py_bases[i].linker;
        if (vec.size() != 2) {
            cerr << "Incorrect number of elements specified for Field \"Backbone_Connect\" "
                 << "\". There must be exactly 2 indices indicated." << endl;
            exit(1);
        }
        array<size_t, 2> linkers = {vec[0], vec[1]};

        auto name = py_bases[i].name;
        transform(name.begin(), name.end(), name.begin(), ::tolower);
        auto code = py_bases[i].code;
        transform(code.begin(), code.end(), code.begin(), ::toupper);
        auto pair_name = py_bases[i].pair_name;
        transform(pair_name.begin(), pair_name.end(), pair_name.begin(), ::tolower);
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
    vector<Base> bases;
    for (auto &v : strand) {
        transform(v.begin(), v.end(), v.begin(), ::tolower);
        bases.push_back(getBaseFromName(v));
    }
    return bases;
}

std::vector<Base> Bases::getComplimentBasesFromStrand(std::vector<std::string> strand) {
    if (all_bases_pair) {
        vector<Base> compliment_bases;
        for (auto &v : strand) {
            transform(v.begin(), v.end(), v.begin(), ::tolower);
            compliment_bases.push_back(name_base_map.find(v)->second);
        }
        return compliment_bases;
    }
    return vector<Base>();
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

    backbone.rotate(rot);
    backbone.translate(atoms[0]->GetVector() - atoms[3]->GetVector());

    // Deleting old hydrogens that formed the vector
    backbone.deleteVectorAtom();
    base.deleteVectorAtom();

    OBMol mol;
    mol += base.getMolecule();
    base_ += base.getMolecule();

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
        conv.Write(&mol);
    }

    perceiveN3OrN1();

    unit = mol;
    base_connect_index = num_atoms + backbone.getLinker()->GetIdx();
    backbone_interconnects = {backbone.getHead()->GetIdx() + num_atoms, backbone.getTail()->GetIdx() + num_atoms};
    base_index_range = {1, num_atoms};
    backbone_index_range = {num_atoms + 1, mol.NumAtoms()};
}

void BaseUnit::perceiveN3OrN1() {
    vector3 zero{0, 0, 0};
    pair<double,unsigned> n1_dist, n3_dist;
    unsigned c6_index;
    auto res = base_.GetResidue(0);
    for (auto it = res->BeginAtoms(); it != res->EndAtoms(); ++it) {
        auto a = *it;
        auto a_num = a->GetAtomicNum();
        auto a_type = string(a->GetType());
        a_type = res->GetAtomID(a);
        transform(a_type.begin(), a_type.end(), a_type.begin(), ::toupper);
        if (a_num == 7 /*nitgrogen*/) {
            if (a_type.find("N3") != string::npos) {
                n3_dist = pair<double,unsigned>(a->GetVector().distSq(zero), a->GetIdx());
            } else if (a_type.find("N1") != string::npos) {
                n1_dist = pair<double,unsigned>(a->GetVector().distSq(zero), a->GetIdx());
            }
        } else if (a_num == 6) {
            if (a_type.find("C6") != string::npos) {
                c6_index = a->GetIdx();
            }
        }
    }
    if (n1_dist.first > n3_dist.first) {
        c6_to_n3 = tuple<string,unsigned,unsigned>("N3", n3_dist.second, c6_index);
    } else if (n3_dist.first > n1_dist.first) {
        c6_to_n3 = tuple<string,unsigned,unsigned>("N1", n1_dist.second, 0);
    } else {
        cerr << "There was an error in perceiving N3 or N1 hydrogen bonding..." << endl;
        throw(1);
    }
}

void RuntimeParameters::validate() {

    if (!std::isnan(base_to_backbone_bond_length) & base_to_backbone_bond_length < 0) {
        cerr << "Base_to_Backbone_Bond_Length must be greater than 0." << endl;
        exit(1);
    }
    if (max_distance < 0) {
        cerr << "Max_Distance must be greater than 0." << endl;
        exit(1);
    }
}

HelicalParameters::HelicalParameters(PyHelicalParameters &hp) {

    vector<RangeField*> ranges{&tilt, &roll, &twist, &shift, &slide, &rise, &buckle, &propeller,
                               &opening, &shear, &stretch, &stagger, &inclination, &tip,
                               &x_displacement, &y_displacement};

    tilt.d = hp.tilt;
    roll.d = hp.roll;
    twist.d = hp.twist;
    shift.d = hp.shift;
    slide.d = hp.slide;
    rise.d = hp.rise;
    buckle.d = hp.buckle;
    propeller.d = hp.propeller;
    opening.d = hp.opening;
    shear.d = hp.shear;
    stretch.d = hp.stretch;
    stagger.d = hp.stagger;
    inclination.d = hp.inclination;
    tip.d = hp.tip;
    x_displacement.d = hp.x_displacement;
    y_displacement.d = hp.y_displacement;

    unsigned i = 0;
    for (auto &r : ranges) {
        checkRangeField(*r);
    }

}
