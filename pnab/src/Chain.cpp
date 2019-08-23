//
// Created by jbarnett8 on 2/10/18.
//

#include <openbabel/obconversion.h>
#include "Chain.h"

using namespace std;
using namespace PNAB;
using namespace OpenBabel;

Chain::Chain(Bases bases, const Backbone &backbone, vector<string> strand, string ff_type,
             std::array<unsigned, 2> &range, bool double_stranded, bool hexad, vector<bool> strand_orientation) {

    strand_orientation_ = strand_orientation;
    hexad_ = hexad;
    double_stranded_ = double_stranded;
    ff_type_ = ff_type;
    monomer_bb_index_range_ = range;
    backbone_ = backbone;
    chain_length_ = static_cast<unsigned>(strand.size());
    pFF_ = OBForceField::FindForceField(ff_type_);
    if (!pFF_) {
        cerr << "Cannot find force field. Exiting" << endl;
        exit(1);
    }
    isKCAL_ = pFF_->GetUnit().find("kcal") != string::npos;

    auto bases_A = bases.getBasesFromStrand(strand);
    auto bases_B = bases.getComplimentBasesFromStrand(strand);

    canonical_nucleobase_ = false;
    for (auto v: bases_A) {
        if (canonical_nucleobase_)
            break;
        for (auto i: std::vector<string> {"a", "g", "c", "t", "u"}) {
            if (v.getName() == i) {
                canonical_nucleobase_ = true;
                break;
            }
        }
    }

    std::vector<std::vector<Base>> v_bases{bases_A, bases_B};
    if (bases_B.empty() && (double_stranded_ || hexad_)) {
        cerr << "User requested double stranded or hexad chain, however not all bases have a valid compliment. Please check"
             << " input file. Exiting..." << endl;
        throw 1;
    }

    n_chains_ = 1;
    if (double_stranded_)
        n_chains_ = 2;
    else if (hexad_)
        n_chains_ = 6;
    unsigned offset = 0;
    for (unsigned i=0; i < n_chains_; i++) {
        setupChain(v_bases[i%2], v_chain_[i], v_new_bond_ids_[i], v_deleted_atoms_ids_[i], v_num_bu_A_mol_atoms_[i], v_bb_start_index_[i],
                   v_base_coords_vec_[i], v_fixed_bonds[i], i);
        std::sort(v_deleted_atoms_ids_[i].begin(), v_deleted_atoms_ids_[i].end());
        setupFFConstraints(v_chain_[i], v_new_bond_ids_[i], v_fixed_bonds[i], offset);
        offset += v_chain_[i].NumAtoms();
        combined_chain_ += v_chain_[i];
    }
}


ConformerData Chain::generateConformerData(double *conf, HelicalParameters &hp, vector<double> energy_filter) {
    if (v_chain_[0].NumAtoms() <= 0) {
        cout << "There are no atoms in the chain_. Exiting..." << endl;
        exit(0);
    }

    unsigned num_cooords;
    if (double_stranded_)
        num_cooords = (v_chain_[0].NumAtoms() + v_chain_[1].NumAtoms()) * 3;
    else if (hexad_)
        num_cooords = (v_chain_[0].NumAtoms() * 3 + v_chain_[1].NumAtoms() * 3) * 3;
    else
        num_cooords = v_chain_[0].NumAtoms() * 3;
    auto *xyz = new double[num_cooords];


    for (unsigned i=0; i < n_chains_; i++) {
        setCoordsForChain(xyz,conf,hp,v_num_bu_A_mol_atoms_[i],v_bb_start_index_[i], v_base_coords_vec_[i],v_deleted_atoms_ids_[i], i);
    }

    ConformerData data;
    data.coords = xyz;
    data.chain_coords_present = true;
    fillConformerEnergyData(xyz, data, energy_filter);
    return data;
}

void Chain::fillConformerEnergyData(double *xyz, PNAB::ConformerData &conf_data, vector<double> energy_filter) {

    OBMol *current_mol;

    current_mol = &combined_chain_;

    current_mol->SetCoordinates(xyz);

    auto n = chain_length_;
    int subtract = 1;
    if (double_stranded_) {
        n *= 2;
        subtract = 2;
    }
    else if (hexad_) {
        n *= 6;
        subtract = 6;
    }

    if (!(n == 1) && !(double_stranded_ && n == 2) && !(hexad_ && n==6)) {

        // Get bond energy
        pFF_->Setup(*current_mol, constraintsBond_);
        conf_data.bondE = pFF_->E_Bond(false) / (n - subtract);

        if (!isKCAL_)
            conf_data.bondE *= KJ_TO_KCAL;

        if (energy_filter[0] < conf_data.bondE) {
            conf_data.accepted = false;
            return;
        }
        // Get angle energy
        pFF_->Setup(*current_mol, constraintsAng_);
        conf_data.angleE = pFF_->E_Angle(false) / (n - subtract);

        if (!isKCAL_)
            conf_data.angleE *= KJ_TO_KCAL;

        if (energy_filter[1] < conf_data.angleE) {
            conf_data.accepted = false;
            return;
        }
    }

    else {
        conf_data.bondE = conf_data.angleE = 0.0;
    }

    // Get torsion energy
    // Add interaction groups for rotatable torsions
    for (int i=0; i < all_torsions_.size(); i++) {
        if (is_fixed_bond[i])
            continue;
        OBBitVec bit = OBBitVec();
        for (auto j: all_torsions_[i]) {
            bit.SetBitOn(j);
        }
        pFF_->AddIntraGroup(bit);
    }

    // We need to trick openbabel to do a new setup for the force field
    // We change the atomic number of one element and then we change it back
    // This triggers a new force field setup
    // This is necessary because for some reason we cannot set interaction
    // groups and then delete them. If no new setup was performed, openbabel
    // would not compute van der Waals or total energy terms for terms not
    // included within the interaction group.
    // We need to keep the torsional computation before van der Waals
    // because it is less expensive. We check for energy thresholds sequentially.
    // Maybe there is a better solution
    unsigned a0 = current_mol->GetAtom(1)->GetAtomicNum();
    current_mol->GetAtom(1)->SetAtomicNum(1);
    pFF_->Setup(*current_mol, constraintsTor_);
    current_mol->GetAtom(1)->SetAtomicNum(a0);
    pFF_->Setup(*current_mol, constraintsTor_);

    conf_data.torsionE = pFF_->E_Torsion(false)/n;
    // This is still necessary
    pFF_->ClearGroups();

    if (!isKCAL_)
        conf_data.torsionE *= KJ_TO_KCAL;

    if (energy_filter[2] < conf_data.torsionE) {
        conf_data.accepted = false;
        return;
    }
    
    // Get torsion energy for fixed bonds
    // Add interaction groups for fixed rotatable torsions
    // Add an empty OBBitVec in case there are not any fixed bonds
    OBBitVec bit = OBBitVec();
    pFF_->AddIntraGroup(bit);
    for (int i=0; i < all_torsions_.size(); i++) {
        if (!is_fixed_bond[i])
            continue;
        OBBitVec bit = OBBitVec();
        for (auto j: all_torsions_[i]) {
            bit.SetBitOn(j);
        }
        pFF_->AddIntraGroup(bit);
    }

    // Do the same trick
    unsigned a1 = current_mol->GetAtom(1)->GetAtomicNum();
    current_mol->GetAtom(1)->SetAtomicNum(1);
    pFF_->Setup(*current_mol, constraintsTor_);
    current_mol->GetAtom(1)->SetAtomicNum(a1);
    pFF_->Setup(*current_mol, constraintsTor_);

    pFF_->E_Torsion(false)/n;
    conf_data.fixed_torsionE = pFF_->E_Torsion(false)/n;
    // This is still necessary
    pFF_->ClearGroups();

    if (!isKCAL_)
       conf_data.fixed_torsionE *= KJ_TO_KCAL;


    // Get VDW energy
    unsigned a = current_mol->GetAtom(1)->GetAtomicNum();
    current_mol->GetAtom(1)->SetAtomicNum(1);
    pFF_->Setup(*current_mol, constraintsTot_);
    current_mol->GetAtom(1)->SetAtomicNum(a);
    pFF_->Setup(*current_mol, constraintsTot_);

    conf_data.VDWE = pFF_->E_VDW(false) / n;

    if (!isKCAL_)
        conf_data.VDWE *= KJ_TO_KCAL;

    if (energy_filter[3] < conf_data.VDWE) {
        conf_data.accepted = false;
        return;
    }

    // Get total energy
    conf_data.total_energy = pFF_->Energy(false) / n;

    if (!isKCAL_)
        conf_data.total_energy *= KJ_TO_KCAL;

    if (energy_filter[4] < conf_data.total_energy) {
        conf_data.accepted = false;
        return;
    }

    conf_data.accepted = true;
}

void Chain::setupChain(std::vector<PNAB::Base> &strand, OpenBabel::OBMol &chain, std::vector<unsigned> &new_bond_ids,
                       std::vector<unsigned> &deleted_atoms_ids, std::vector<unsigned> &num_base_unit_atoms,
                       std::vector<unsigned> &bb_start_index, std::vector<double *> &base_coords_vec,
                       std::vector<std::vector<unsigned>> &fixed_bonds_vec, unsigned chain_index) {

    auto bases_A = strand;
    unsigned neighbor_id;
    vector<unsigned> head_ids, tail_ids;
    OBAtom *head, *tail;

    vector<BaseUnit> bu_A;
    vector<string> base_codes;
    vector<string> base_names;
    for (auto v : bases_A) {
        bu_A.emplace_back(BaseUnit(v, backbone_));
        base_codes.emplace_back(v.getCode());
        base_names.emplace_back(v.getName());
    }

    vector<OBMol> bu_A_mol;
    vector<array<size_t, 2>> bu_linkers_vec;
    vector<vector<vector<unsigned>>> base_fixed_bonds;
    for (auto v : bu_A) {
        // We want a list of linkers for deletion
        bu_linkers_vec.emplace_back(v.getBackboneLinkers());

        // Place where
        //base_connect_index.emplace_back(v.getBaseConnectIndex());
        // We need to know when the backbone starts in the BaseUnit
        auto bb_r = v.getBackboneIndexRange();
        bb_start_index.emplace_back(static_cast<unsigned>(bb_r[0]));

        // We want some information from each of the molecules, specifically the size
        auto mol = v.getMol();
        num_base_unit_atoms.emplace_back(mol.NumAtoms());
        bu_A_mol.emplace_back(mol);
        auto *xyz = new double[num_base_unit_atoms.back() * 3];
        memcpy(xyz,mol.GetCoordinates(),sizeof(double) * num_base_unit_atoms.back() * 3);
        base_coords_vec.emplace_back(xyz);

        // Fixed bonds
        base_fixed_bonds.push_back(v.getFixedBonds());
    }

    auto chain_letter = static_cast<char>('A' + chain_index);
    unsigned c = 0;
    for (auto v : bu_A_mol) {
        for (auto it = v.BeginResidues(); it != v.EndResidues(); ++it) {
            OBResidue *r = *it;
            r->SetName(base_codes[c]);
            r->SetChainNum(chain_index);
            r->SetChain(chain_letter);
            r->SetTitle(base_names[c].c_str());
            r->SetNum(c + 1);
        }
        chain += v;
        c++;
    }

    // Mark fixed bonds
    int num_atoms = 0, n_fixed = 0;
    for (int i=0; i < base_fixed_bonds.size(); i++) {
        for (auto bonds: base_fixed_bonds[i]) {
            OBAtom* a1 = chain.GetAtom(bonds[0] + num_atoms);
            OBAtom* a2 = chain.GetAtom(bonds[1] + num_atoms);
            OBPairData *label1 = new OBPairData();
            OBPairData *label2 = new OBPairData();
            label1->SetAttribute(to_string(n_fixed));
            label2->SetAttribute(to_string(n_fixed));
            a1->SetData(label1);
            a2->SetData(label2);
            n_fixed++;
        }
        num_atoms += num_base_unit_atoms[i];
    }

    bu_linkers_vec.pop_back();

    std::vector<long int> center_ids;
    unsigned count = 0, index = 0;
    if ((canonical_nucleobase_ && chain_index == 0) || (!canonical_nucleobase_ && strand_orientation_[chain_index])) {
        for (unsigned i = 0; i < bu_linkers_vec.size(); ++i) {
            head_ids.push_back(static_cast<int>(chain.GetAtom(count + bu_linkers_vec[i][0])->GetId()));
            count += num_base_unit_atoms[index++];
            tail_ids.push_back(static_cast<int>(chain.GetAtom(count + bu_linkers_vec[i + 1][1])->GetId()));
        }
    }

    else {
        for (unsigned i = 0; i < bu_linkers_vec.size(); ++i) {
            tail_ids.push_back(static_cast<int>(chain.GetAtom(count + bu_linkers_vec[i    ][1])->GetId()));
            count += num_base_unit_atoms[index++];
            head_ids.push_back(static_cast<int>(chain.GetAtom(count + bu_linkers_vec[i + 1][0])->GetId()));
        }
    }

    for (unsigned i = 0; i < head_ids.size(); i++) {
        head = chain.GetAtomById(head_ids[i]);
        tail = chain.GetAtomById(tail_ids[i]);
        FOR_NBORS_OF_ATOM(nbr, tail) {
            if (nbr->GetAtomicNum() == 1) {
                deleted_atoms_ids.push_back(nbr->GetId());
                chain.DeleteAtom(chain.GetAtom(nbr->GetIdx()));
                break;
            }
        }
        FOR_NBORS_OF_ATOM(nbr, head) {
            if (nbr->GetAtomicNum() == 1) {
                deleted_atoms_ids.push_back(nbr->GetId());
                chain.DeleteAtom(chain.GetAtom(nbr->GetIdx()));
                break;
            }
        }
        FOR_NBORS_OF_ATOM(nbr, tail) {
            neighbor_id = nbr->GetId();
            break;
        }
        deleted_atoms_ids.push_back(tail->GetId());
        chain.DeleteAtom(tail);
        new_bond_ids.push_back(head->GetId());
        new_bond_ids.push_back(neighbor_id);
        chain.AddBond(head->GetIdx(), chain.GetAtomById(neighbor_id)->GetIdx(), 1);
    }


    for (int i=0; i < n_fixed; i++) {
        fixed_bonds_vec.push_back({});
        FOR_ATOMS_OF_MOL(a, chain) {
            if (a->HasData(to_string(i))) {
                fixed_bonds_vec[i].push_back(a->GetIdx());
            }
        }
    }

};

void Chain::setupFFConstraints(OpenBabel::OBMol &chain, std::vector<unsigned> &new_bond_ids, std::vector<std::vector<unsigned>> &fixed_bonds_vec, unsigned offset) {
    int nbrIdx;

    // Torsional energy, must fix non-torsion atoms
    vector<int> torsionAtoms;
    for (auto v : new_bond_ids) {
        FOR_NBORS_OF_ATOM(nbr, chain.GetAtomById(v))
            torsionAtoms.emplace_back(nbr->GetIdx());
    }
    FOR_ATOMS_OF_MOL(a, chain) {
        if(std::find(torsionAtoms.begin(), torsionAtoms.end(), a->GetIdx()) == torsionAtoms.end()) {
            constraintsAng_.AddIgnore(a->GetIdx() + offset);
            constraintsBond_.AddIgnore(a->GetIdx() + offset);
        }
    }

    // Angle energy, must fix non-angle atoms
    for (int i = 1; i < new_bond_ids.size(); i+=2) {
        FOR_NBORS_OF_ATOM(nbr, chain.GetAtomById(new_bond_ids.at(i))) {
            nbrIdx = nbr->GetIdx();
            if (nbrIdx != chain.GetAtomById(new_bond_ids.at(i - 1))->GetIdx()) {
                constraintsAng_.AddIgnore(nbrIdx + offset);
                constraintsBond_.AddIgnore(nbrIdx + offset);
            }
        }
    }

    // Bond energy, must fix non-bond atoms
    for (int i = 0; i < new_bond_ids.size(); i+=2) {
        FOR_NBORS_OF_ATOM(nbr, chain.GetAtomById(new_bond_ids.at(i))) {
            nbrIdx = nbr->GetIdx();
            if (nbrIdx != chain.GetAtomById(new_bond_ids.at(i + 1))->GetIdx()) {
                constraintsBond_.AddIgnore(nbrIdx + offset);
            }
        }
    }

    // Determine torsion atoms
    OBRotorList rl;
    OBBitVec fix_bonds(chain.NumAtoms());
    rl.Setup(chain);
    rl.SetRotAtomsByFix(chain);
    OBRotorIterator ri;
    OBRotor *r = rl.BeginRotor(ri);

    while(r) {
        vector<int> v = r->GetDihedralAtoms();
        FOR_NBORS_OF_ATOM(nbr, chain.GetAtom(v[1])) {
            if (nbr->GetIdx() == v[2])
                continue;
            FOR_NBORS_OF_ATOM(nbr2, chain.GetAtom(v[2])) {
                if (nbr2->GetIdx() == v[1])
                    continue;
                unsigned num = nbr->GetResidue()->GetNum();
                // Exclude dihedral angles in two different residues as these bonds are not
                // rotated in the conformtation search 
                if ((chain.GetAtom(v[1])->GetResidue()->GetNum() == num) &&
                    (chain.GetAtom(v[2])->GetResidue()->GetNum() == num) &&
                    (nbr2->GetResidue()->GetNum() == num)) {
                    // These are the correct torsional angles that we want to compute energies for
                    bool fixed = false;
                    for (auto bond: fixed_bonds_vec) {
                        if ((find(bond.begin(), bond.end(), v[1]) != bond.end()) && (find(bond.begin(), bond.end(), v[2]) != bond.end())) {
                            is_fixed_bond.push_back(true);
                            fixed = true;
                            break;
                        }
                    }
                    if (!fixed)
                        is_fixed_bond.push_back(false);
                    all_torsions_.push_back({nbr->GetIdx() + offset, v[1] + offset, v[2] + offset, nbr2->GetIdx() + offset});
                }
            }
        }

        r = rl.NextRotor(ri);
    }
}

void Chain::setCoordsForChain(double *xyz, double *conf, PNAB::HelicalParameters &hp,
                              std::vector<unsigned> &num_bu_atoms,
                              std::vector<unsigned> &bb_start_index, std::vector<double *> &base_coords_vec,
                              std::vector<unsigned> &deleted_atoms_ids, unsigned chain_index) {

    matrix3x3 change_sign, z_rot;
    matrix3x3 change_sign2 = {{-1, 0, 0},{0, 1, 0}, {0, 0, -1}};
    auto g_rot = hp.getGlobalRotationOBMatrix();
    auto g_trans = hp.getGlobalTranslationVec();
    unsigned xyzI = 0, local_offset = 0, deleted_atom_index = 0;
    if (double_stranded_ && chain_index == 1) {
        xyzI = 3 * v_chain_[0].NumAtoms();
        change_sign = {{1, 0, 0},{0, -1, 0}, {0, 0, -1}};
    }
    else if (!canonical_nucleobase_) {
        if (hexad_)
            xyzI = 3 * v_chain_[0].NumAtoms() * ((chain_index + 1) / 2) + 3 * v_chain_[1].NumAtoms() * (chain_index / 2); 
        change_sign = {{1, 0, 0},{0, -1, 0}, {0, 0, -1}};
        double twist = chain_index * 60.0 * DEG_TO_RAD;
        z_rot = {{cos(twist), -sin(twist), 0},{sin(twist), cos(twist), 0}, {0, 0, 1}};
    }

    for (unsigned i = 0; i < chain_length_; ++i) {
        auto n = num_bu_atoms[i];
        auto r = bb_start_index[i];
        auto s_trans = hp.getStepTranslationVec(i);
        auto s_rot   = hp.getStepRotationOBMatrix(i);
        double *base_coords = base_coords_vec[i];
        for (unsigned baseI = 0; baseI < r - 1; ++baseI) {
            vector3 v3;
            v3.Set(base_coords + 3 * baseI);

            v3 += g_trans; v3 *= g_rot;

            if (double_stranded_ && chain_index == 1)
                v3 *= change_sign;

            else if (!canonical_nucleobase_) {
                if (!strand_orientation_[chain_index])
                    v3 *= change_sign;
                v3 *= z_rot;
            }


            v3 += s_trans; v3 *= s_rot;

            // Reflect to get the orientation that agrees with 3DNA for DNA/RNA
            //v3 *=  change_sign2;

            v3.Get(xyz + xyzI);
            xyzI += 3;
        }
        unsigned monomer_index = monomer_bb_index_range_[0] - 1;
        for (unsigned bbI = r - 1; bbI < n; bbI++) {
            if (deleted_atoms_ids.empty() || bbI + local_offset + 1 != deleted_atoms_ids[deleted_atom_index] + 1) {
                vector3 v3;
                v3.Set(conf + 3 * monomer_index);

                v3 += g_trans; v3 *= g_rot;

                if (double_stranded_ && chain_index == 1)
                    v3 *= change_sign;

                else if (!canonical_nucleobase_) {
                    if (!strand_orientation_[chain_index])
                        v3 *= change_sign;
                    v3 *= z_rot;
                }
               
                v3 += s_trans; v3 *= s_rot;

                // Reflect to get the orientation that agrees with 3DNA for DNA/RNA
                //v3 *=  change_sign2;

                v3.Get(xyz + xyzI);
                xyzI += 3;
            } else {
                deleted_atom_index++;
            }
            monomer_index++;
        }
        local_offset += n;
    }

}
