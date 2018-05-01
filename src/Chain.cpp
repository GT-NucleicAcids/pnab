//
// Created by jbarnett8 on 2/10/18.
//

#include <openbabel/obconversion.h>
#include "Chain.h"

using namespace std;
using namespace PNAB;
using namespace OpenBabel;

Chain::Chain(Bases bases, Backbone backbone, vector<string> strand, string ff_type,
             std::array<unsigned, 2> &range, bool double_stranded) {

    unsigned neighborId;
    vector<unsigned> headIds, tailIds;
    OBAtom *head, *tail;

    // New code -----------------------------------------------------------------------------------
    monomer_bb_index_range_ = range;
    auto bases_A = bases.getBasesFromStrand(strand);
    chain_length_ = static_cast<unsigned>(strand.size());
    vector<BaseUnit> bu_A;
    for (auto v : bases_A)
        bu_A.push_back(BaseUnit(v,backbone));
    vector<OBMol> bu_A_mol;
    vector<array<size_t, 2>> bu_linkers_vec;
    for (auto v : bu_A) {
        // We want a list of linkers for deletion
        bu_linkers_vec.push_back(v.getBackboneLinkers());

        // We need to know when the backbone starts in the BaseUnit
        auto bb_r = v.getBackboneIndexRange();
        bb_start_index_.push_back(static_cast<unsigned>(bb_r[0]));

        // We want some information from each of the molecules, specifically the size
        auto mol = v.getMol();
        num_bu_A_mol_atoms_.push_back(mol.NumAtoms());
        bu_A_mol.push_back(mol);
        double *xyz = new double[num_bu_A_mol_atoms_.back() * 3];
        memcpy(xyz,mol.GetCoordinates(),sizeof(double) * num_bu_A_mol_atoms_.back() * 3);
        base_coords_vec_.push_back(xyz);
    }

    vector3 v3(0,0,10);
    unsigned c = 0;
    for (auto v : bu_A_mol) {
        v.Translate(v3 * c);
        chain_ += v;
        c++;
    }
    

    bu_linkers_vec.pop_back();

    unsigned count = 0, index = 0;
    for (unsigned i = 0; i < bu_linkers_vec.size(); ++i) {
        headIds.push_back(static_cast<int>(chain_.GetAtom(count + bu_linkers_vec[i    ][0])->GetId()));
        count += num_bu_A_mol_atoms_[index++];
        tailIds.push_back(static_cast<int>(chain_.GetAtom(count + bu_linkers_vec[i + 1][1])->GetId()));
    }

//    cout << "Heads and tails---" << endl;
//    for (auto v : headIds)
//        cout << v + 1<< ", ";
//    cout << endl;
//    for (auto v : tailIds)
//        cout << v + 1<< ", ";
//    cout << endl;

    for (unsigned i = 0; i < headIds.size(); i++) {
        head = chain_.GetAtomById(headIds[i]);
        tail = chain_.GetAtomById(tailIds[i]);
        FOR_NBORS_OF_ATOM(nbr, tail) {
            if (nbr->IsHydrogen()) {
                deletedAtomsId_.push_back(nbr->GetId());
                chain_.DeleteAtom(chain_.GetAtom(nbr->GetIdx()));
                break;
            }
        }
        FOR_NBORS_OF_ATOM(nbr, head) {
            if (nbr->IsHydrogen()) {
                deletedAtomsId_.push_back(nbr->GetId());
                chain_.DeleteAtom(chain_.GetAtom(nbr->GetIdx()));
                break;
            }
        }
        FOR_NBORS_OF_ATOM(nbr, tail) {
            neighborId = nbr->GetId();
            break;
        }
        deletedAtomsId_.push_back(tail->GetId());
        chain_.DeleteAtom(tail);
        newBondIDs_.push_back(head->GetId()); newBondIDs_.push_back(neighborId);
        chain_.AddBond(head->GetIdx(), chain_.GetAtomById(neighborId)->GetIdx(), 1);
    }

//    BaseUnit bu(bases.getBaseFromName(strand[0]),backbone);
//    monomer_ = bu.getMol();
//    unsigned numMonomerAtoms = monomer_.NumAtoms();
//    chain_length_ = strand.size();
//
//    chain_ = OBMol();
//    auto bu_linkers = bu.getBackboneLinkers();
//
//    // Creating raw chain_
//    for (unsigned polymerI = 0; polymerI < chain_length_; polymerI++) {
//        chain_+=monomer_;
//    }
//
//    // Deleting appropriate atoms and creating bonds
//    for (int polymerI = 0; polymerI < chain_length_ - 1; polymerI++) {
//        headIds.push_back(chain_.GetAtom(numMonomerAtoms *   polymerI   + bu_linkers[0])->GetId());
//        tailIds.push_back(chain_.GetAtom(numMonomerAtoms * (polymerI+1) + bu_linkers[1])->GetId());
//    }


//    for (int polymerI = 0; polymerI < chain_length_ - 1; polymerI++) {
//        head = chain_.GetAtomById(headIds.at(polymerI));
//        tail = chain_.GetAtomById(tailIds.at(polymerI));
//        FOR_NBORS_OF_ATOM(nbr, tail) {
//            if (nbr->IsHydrogen()) {
//                deletedAtomsId_.push_back(nbr->GetId());
//                chain_.DeleteAtom(chain_.GetAtom(nbr->GetIdx()));
//                break;
//            }
//        }
//        FOR_NBORS_OF_ATOM(nbr, head) {
//            if (nbr->IsHydrogen()) {
//                deletedAtomsId_.push_back(nbr->GetId());
//                chain_.DeleteAtom(chain_.GetAtom(nbr->GetIdx()));
//                break;
//            }
//        }
//        // chain_.DeleteHydrogens(tail);  chain_.DeleteHydrogens(head);
//        FOR_NBORS_OF_ATOM(nbr, tail) {
//            neighborId = nbr->GetId();
//            break;
//        }
//        deletedAtomsId_.push_back(tail->GetId());
//        chain_.DeleteAtom(tail);
//        newBondIDs_.push_back(head->GetId()); newBondIDs_.push_back(neighborId);
//        chain_.AddBond(head->GetIdx(), chain_.GetAtomById(neighborId)->GetIdx(), 1);
//    }

    std::sort (deletedAtomsId_.begin(), deletedAtomsId_.end());

    pFF_ = OBForceField::FindForceField("GAFF");
    pFF_->Setup(chain_);
    isKCAL_ = pFF_->GetUnit().find("kcal") != string::npos;

    // Setting up constraints
    constraintsTot_ = new OBFFConstraints;
    constraintsTor_ = new OBFFConstraints;
    constraintsAng_ = new OBFFConstraints;
    constraintsBond_ = new OBFFConstraints;
    int nbrIdx;

    // Torsional energy, must fix non-torsion atoms
    vector<int> torsionAtoms;
    for (int i = 0; i < newBondIDs_.size(); i++) {
        FOR_NBORS_OF_ATOM(nbr, chain_.GetAtomById(newBondIDs_.at(i))) {
            torsionAtoms.push_back(nbr->GetIdx());
        }
    }
    FOR_ATOMS_OF_MOL(a, chain_) {
        if(std::find(torsionAtoms.begin(), torsionAtoms.end(), a->GetIdx()) == torsionAtoms.end()) {
            constraintsTor_->AddIgnore(a->GetIdx());
            constraintsAng_->AddIgnore(a->GetIdx());
            constraintsBond_->AddIgnore(a->GetIdx());
        }
    }

    // Angle energy, must fix non-angle atoms
    for (int i = 1; i < newBondIDs_.size(); i+=2) {
        FOR_NBORS_OF_ATOM(nbr, chain_.GetAtomById(newBondIDs_.at(i))) {
            nbrIdx = nbr->GetIdx();
            if (nbrIdx != chain_.GetAtomById(newBondIDs_.at(i - 1))->GetIdx()) {
                constraintsAng_->AddIgnore(nbrIdx);
                constraintsBond_->AddIgnore(nbrIdx);
            }
        }
    }

    // Bond energy, must fix non-bond atoms
    for (int i = 0; i < newBondIDs_.size(); i+=2) {
        FOR_NBORS_OF_ATOM(nbr, chain_.GetAtomById(newBondIDs_.at(i))) {
            nbrIdx = nbr->GetIdx();
            if (nbrIdx != chain_.GetAtomById(newBondIDs_.at(i + 1))->GetIdx()) {
                constraintsBond_->AddIgnore(nbrIdx);
            }
        }
    }

    if (!(constraintsTot_ && constraintsTor_ && constraintsAng_ && constraintsBond_)) {
        cerr << "Error initializing force field constraints for chain_. Exiting..." << endl;
        exit(1);
    }


    // garbage code for debugging
    // TODO delete this code
    {
        OBConversion conv_;
        std::filebuf fb;
        fb.open("test_chain_connected.pdb", std::ios::out);
        std::ostream fileStream(&fb);
        conv_.SetOutFormat("PDB");
        conv_.SetOutStream(&fileStream);
        conv_.Write(&chain_);
    }

}

ConformerData Chain::generateConformerData(double *conf, HelicalParameters &hp) {
    if (chain_.NumAtoms() <= 0) {
        cout << "There are no atoms in the chain_. Exiting..." << endl;
        exit(0);
    }

    // new new *new* code
    auto g_rot = hp.getGlobalRotationOBMatrix();
    auto g_trans = hp.getGlobalTranslationVec();
    double *xyz = new double[chain_.NumAtoms() * 3];
    unsigned xyzI = 0, offset = 0, deleted_atom_index = 0;

//    cout << "DELETED ATOMS::::::::" << endl;
//    for (auto v : deletedAtomsId_)
//        cout << v + 1 << ", ";
//    cout << endl << endl;

    for (unsigned i = 0; i < chain_length_; ++i) {
        auto n = num_bu_A_mol_atoms_[i];
        auto r = bb_start_index_[i];
        auto s_trans = hp.getStepTranslationVec(i);
        auto s_rot   = hp.getStepRotationOBMatrix(i);
        double *base_coords = base_coords_vec_[i];
        for (unsigned baseI = 0; baseI < r - 1; ++baseI) {
//            cout << "Base: " << baseI + offset + 1 << endl;
            vector3 v3;
            v3.Set(base_coords + 3 * baseI);
//            cout << "Base: baseI->" << baseI << ", " << v3 << endl;
            v3 += g_trans;
            v3 *= g_rot;
            v3 += s_trans;
            v3 *= s_rot;
            v3.Get(xyz + xyzI);
            xyzI += 3;
        }
        unsigned monomer_index = monomer_bb_index_range_[0] - 1;
        for (unsigned bbI = r - 1; bbI < n; bbI++) {
            if (bbI + offset + 1 != deletedAtomsId_[deleted_atom_index] + 1) {
//                cout << "Backbone: " << bbI + offset + 1 << endl;
                vector3 v3;
                v3.Set(conf + 3 * monomer_index);
//                cout << "Backbone: bbI->" << bbI << ", " << v3 << endl;
                v3 += s_trans;
                v3 *= s_rot;
                v3.Get(xyz + xyzI);
                xyzI += 3;
            } else {
//                cout << "Skipped...." << endl;
                deleted_atom_index++;
            }
            monomer_index++;
        }
        offset += n;
    }


    // new new code
//    auto g_rot = hp.getGlobalRotationOBMatrix();
//    auto g_trans = hp.getGlobalTranslationVec();
//
//    double *xyz = new double[chain_.NumAtoms() * 3];
//
//    unsigned xyzI = 0, chain_count = 0, offset = 0, ai = 0;
//    int deleted_atom_index = 0;
//    for (auto n : num_bu_A_mol_atoms_) {
//        auto s_trans = hp.getStepTranslationVec(chain_count);
//        auto s_rot   = hp.getStepRotationOBMatrix(chain_count);
//        for (unsigned i = 0; i < n; ++i) {
//            if (i + offset + 1 != deletedAtomsId_[deleted_atom_index] + 1) {
//                vector3 v;
//                v.Set(conf + 3 * i);
//                v += s_trans;
//                v *= s_rot;
//                v.Get(xyz + xyzI);
//                xyzI += 3;
//            } else {
//                deleted_atom_index++;
//            }
//        }
//        chain_count++;
//        offset += n;
//    }



    // new code
//    auto g_rot = hp.getGlobalRotationOBMatrix();
//    auto g_trans = hp.getGlobalTranslationVec();
//
//    double *xyz = new double[chain_.NumAtoms() * 3];
//
//    unsigned xyzI = 0, chain_count = 0, offset = 0, ai = 0;
//    int deleted_atom_index = 0;
//    for (auto n : num_bu_A_mol_atoms_) {
//        for (unsigned i = 0; i < 3 * n; i += 3) {
//            ai = (i + 1) / 3 + offset + 1;
//            if (deleted_atom_index == -1 || ai != deletedAtomsId_[deleted_atom_index] + 1) {
//                vector3 v;
//                v.Set(conf + i);
//                v += hp.getStepTranslationVec(chain_count);
//                v *= hp.getStepRotationOBMatrix(chain_count);
//                v.Get(xyz + xyzI);
//                xyzI += 3;
////                cout << v << ", " << ai << ", " << deleted_atom_index << endl;
//            } else {
//                if ((deleted_atom_index + 1) < deletedAtomsId_.size())
//                    deleted_atom_index++;
//                else
//                    deleted_atom_index = -1;
//            }
//        }
//        chain_count++;
//        offset += n;
//    }


//    int num_monomer_atoms = monomer_.NumAtoms(), ai, offset, dai = 0, xyzI = -3;
//    double *xyz;
//    auto g_rot = hp.getGlobalRotationOBMatrix();
//    auto g_trans = hp.getGlobalTranslationVec();
//
//    xyz = new double[chain_.NumAtoms() * 3];
//
//    for (int chainIdx = 0; chainIdx < chain_length_; chainIdx++) {
//        for (int ci = 0; ci < 3*num_monomer_atoms; ci+=3) {
//            offset = chainIdx * num_monomer_atoms;
//            ai = (ci + 1) / 3 + offset + 1;
//            /* for some reason, OpenBabel starts id's indexation at 0 instead of 1 like *
//            *  the normal indexation for the atom index. This is the source of the +1  */
//            if (dai == -1 || ai != deletedAtomsId_.at(dai) + 1) {
//                xyzI+=3;
//                vector3 v;
//                v.Set(conf + ci);
//                v += hp.getStepTranslationVec(chainIdx);
//                v *= hp.getStepRotationOBMatrix(chainIdx);
//                v.Get(xyz + xyzI);
////                cout << v << ", " << ai << ", " << dai << endl;
//            } else {
//                if ((dai + 1) < deletedAtomsId_.size())
//                    dai++;
//                else
//                    dai = -1;
//            }
//        }
//    }

    ConformerData data;
    data.coords = xyz;
    data.chain_coords_present = true;
    fillConformerEnergyData(xyz, data);
    return data;
}

void Chain::fillConformerEnergyData(double *xyz, PNAB::ConformerData &conf_data) {
    chain_.SetCoordinates(xyz);

    auto n = chain_length_;

    // Get total energies and VDW
    pFF_->Setup(chain_, *constraintsTot_);
    conf_data.total_energy = pFF_->Energy() / n;
    conf_data.VDWE = pFF_->E_VDW() / n;
    conf_data.totTorsionE = pFF_->E_Torsion() / n;

    if (n - 1 <= 0)
        n = 2;

    // Get torsional energies
    pFF_->SetConstraints(*constraintsTor_);
    conf_data.torsionE = pFF_->E_Torsion() / (n - 1);

    // Get angle energy
    pFF_->SetConstraints(*constraintsAng_);
    conf_data.angleE = pFF_->E_Angle() / (n - 1);

    // Get bond energy
    pFF_->SetConstraints(*constraintsBond_);
    conf_data.bondE = pFF_->E_Bond() / (n - 1);

    if (!isKCAL_) {
        conf_data.total_energy *= KJ_TO_KCAL;
        conf_data.VDWE         *= KJ_TO_KCAL;
        conf_data.totTorsionE  *= KJ_TO_KCAL;
        conf_data.torsionE     *= KJ_TO_KCAL;
        conf_data.angleE       *= KJ_TO_KCAL;
        conf_data.bondE        *= KJ_TO_KCAL;
    }
}
