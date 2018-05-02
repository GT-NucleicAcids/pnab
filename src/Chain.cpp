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

    ff_type_ = ff_type;
    monomer_bb_index_range_ = range;
    auto bases_A = bases.getBasesFromStrand(strand);
    chain_length_ = static_cast<unsigned>(strand.size());
    vector<BaseUnit> bu_A;
    vector<string> base_codes;
    vector<string> base_names;
    for (auto v : bases_A) {
        bu_A.push_back(BaseUnit(v, backbone));
        base_codes.push_back(v.getCode());
        base_names.push_back(v.getName());
    }
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
        for (auto it = v.BeginResidues(); it != v.EndResidues(); ++it) {
            OBResidue *r = *it;
            r->SetName(base_codes[c]);
            r->SetChainNum(1);
            r->SetChain('A');
            r->SetTitle(base_names[c].c_str());
            r->SetNum(c);
        }
//        v.Translate(v3 * c);
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

    std::sort (deletedAtomsId_.begin(), deletedAtomsId_.end());

    pFF_ = OBForceField::FindForceField(ff_type_);
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

    for (unsigned i = 0; i < chain_length_; ++i) {
        auto n = num_bu_A_mol_atoms_[i];
        auto r = bb_start_index_[i];
        auto s_trans = hp.getStepTranslationVec(i);
        auto s_rot   = hp.getStepRotationOBMatrix(i);
        double *base_coords = base_coords_vec_[i];
        for (unsigned baseI = 0; baseI < r - 1; ++baseI) {
            vector3 v3;
            v3.Set(base_coords + 3 * baseI);
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
                vector3 v3;
                v3.Set(conf + 3 * monomer_index);
                v3 += s_trans;
                v3 *= s_rot;
                v3.Get(xyz + xyzI);
                xyzI += 3;
            } else {
                deleted_atom_index++;
            }
            monomer_index++;
        }
        offset += n;
    }

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
