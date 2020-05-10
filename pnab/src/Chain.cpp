/**@file
 * @brief A file for defining various methods for building and evaluating nucleic acid strands
 */

#include "Chain.h"

using namespace std;
using namespace PNAB;
using namespace OpenBabel;

Chain::Chain(Bases bases, const Backbone &backbone, std::vector<std::string> strand, std::string ff_type,
             std::array<unsigned, 2> &range, bool hexad, std::vector<bool> build_strand, std::vector<bool> strand_orientation,
             double glycosidic_bond_distance) {

    // Set some variables
    glycosidic_bond_distance_ = glycosidic_bond_distance;
    build_strand_ = build_strand;
    strand_orientation_ = strand_orientation;
    hexad_ = hexad;
    ff_type_ = ff_type;
    monomer_bb_index_range_ = range;
    chain_length_ = static_cast<unsigned>(strand.size());

    // Get the openbabel force field
    pFF_ = OBForceField::FindForceField(ff_type_);
    // Make sure we have a valid pointer
    if (!pFF_) {
        throw std::runtime_error("Cannot find force field.");
    }
    // Check whether the force field uses kcal/mol or kJ/mol
    isKCAL_ = pFF_->GetUnit().find("kcal") != string::npos;

    // Get the bases for the strand and the complimentary strand for duplexed and hexads
    auto bases_A = bases.getBasesFromStrand(strand);
    auto bases_B = bases.getComplimentBasesFromStrand(strand);
    std::vector<std::vector<Base>> v_bases{bases_A, bases_B};

    n_chains_ = 0;
    for (auto i: build_strand_) {
        if (i)
            n_chains_ += 1;
    }

//    // Determine the number of chains
    unsigned offset = 0;

    // Setup all chains
    for (unsigned i=0; i < 6; i++) {
        if (build_strand_[i]) {
            setupChain(v_bases[i%2], v_chain_[i], v_new_bond_ids_[i], v_deleted_atoms_ids_[i], v_num_bu_A_mol_atoms_[i], v_bb_start_index_[i],
                       v_base_coords_vec_[i], v_fixed_bonds[i], backbone, i);
            std::sort(v_deleted_atoms_ids_[i].begin(), v_deleted_atoms_ids_[i].end());
            setupFFConstraints(v_chain_[i], v_new_bond_ids_[i], v_fixed_bonds[i], offset);
            offset += v_chain_[i].NumAtoms();
            combined_chain_ += v_chain_[i];
        }
    }
}


ConformerData Chain::generateConformerData(double *conf, HelicalParameters &hp, vector<double> energy_filter) {

    // Set the correct number of coordinates
    unsigned num_cooords = 0;
    for (int i=0; i<6; i++) {
        if (build_strand_[i])
            num_cooords += v_chain_[i].NumAtoms() * 3;
    }

    auto *xyz = new double[num_cooords];

    // Set the coordinates
    for (unsigned i=0; i < 6; i++) {
        if (build_strand_[i])
            setCoordsForChain(xyz,conf,hp,v_num_bu_A_mol_atoms_[i],v_bb_start_index_[i], v_base_coords_vec_[i],v_deleted_atoms_ids_[i], i);
    }

    // Fill in the energy data
    ConformerData data;
    fillConformerEnergyData(xyz, data, energy_filter);

    if (data.accepted) {
        data.molecule = combined_chain_;
        data.molecule.SetChainsPerceived();
        // Reorder atoms
        orderResidues(&data.molecule);
    }

    delete[] xyz;

    return data;
}

void Chain::orderResidues(OBMol* molecule) {
    // Reorder residues in the chain; improves the PDB formatting of systems with inverted chains
    // Basically, this function inverts the order of the residues for inverted chains
    vector<int> order;
    int index = 0;

    for (int i=0; i<6; i++) {
        if (!build_strand_[i])
            continue;

        if (strand_orientation_[i]) {
            for (int j=1; j < v_chain_[i].NumAtoms() + 1; j++)
                order.push_back(j + index);
            index += v_chain_[i].NumAtoms();
        }

        else {
            int residue_per_nucleotide = v_chain_[i].NumResidues() / chain_length_;
            for (int j = v_chain_[i].NumResidues() - 1; j > -1; j = j-residue_per_nucleotide) {
                for (int k=residue_per_nucleotide-1; k > -1; k--) {
                    OBResidue* res = v_chain_[i].GetResidue(j-k);
                    FOR_ATOMS_OF_RESIDUE(atom, res)
                        order.push_back(atom->GetIdx() + index);
                }
            }
            index += v_chain_[i].NumAtoms();
        }
    }

    molecule->RenumberAtoms(order);
}

void Chain::fillConformerEnergyData(double *xyz, PNAB::ConformerData &conf_data, vector<double> energy_filter) {

    // Set the molecule
    OBMol *current_mol;
    current_mol = &combined_chain_;
    current_mol->SetCoordinates(xyz);

    // Determine the number of nucleotides in the system
    unsigned n = 0;
    for (auto i: build_strand_) {
        if(i)
            n += chain_length_; 
    }

    // Get bond energy. We set the constraint to ignore all bonds except the 
    // one bond between the first and second nucleotides in the first strand
    // The other bonds will have the same energy
    pFF_->Setup(*current_mol, constraintsBond_);
    conf_data.bondE = pFF_->E_Bond(false);

    if (!isKCAL_)
        conf_data.bondE *= KJ_TO_KCAL;

    // if not accepted, return
    if (energy_filter[0] < conf_data.bondE) {
        conf_data.accepted = false;
        return;
    }

    // Get angle energy. We compute the energy only between the first and 
    // second nucleotides in the first strand. The other angles will have the same energy

    // We use energy groups here
    // Set the energy groups for the required angles
    OBBitVec bit = OBBitVec();
    pFF_->AddIntraGroup(bit); // Add empty bit in case only one nucleotide per strand is requested
    for (auto i: all_angles_) {
        OBBitVec bit = OBBitVec();
        for (auto j: i)
            bit.SetBitOn(j);
        pFF_->AddIntraGroup(bit);
    }


    pFF_->Setup(*current_mol, constraintsAng_); // empty constraint
    conf_data.angleE = pFF_->E_Angle(false);

    pFF_->ClearGroups();

    if (!isKCAL_)
        conf_data.angleE *= KJ_TO_KCAL;

    // if not accepted, return
    if (energy_filter[1] < conf_data.angleE) {
        conf_data.accepted = false;
        return;
    }

    // Get torsion energy
    // Add interaction groups for rotatable torsions
    bool is_there_fixed_bond = false;
    for (int i=0; i < all_torsions_.size(); i++) {
        // Determine whether we have fixed bonds
        if (is_fixed_bond[i]) {
            is_there_fixed_bond = true;
            continue;
            }
        OBBitVec bit = OBBitVec();
        for (auto j: all_torsions_[i]) {
            bit.SetBitOn(j);
        }
        pFF_->AddIntraGroup(bit);
    }

    constraintsTor_.AddIgnore(0); // A trick to trigger a setup with the new interaction groups
    pFF_->Setup(*current_mol, constraintsTor_);

    conf_data.torsionE = pFF_->E_Torsion(false)/n; // Divide by the number of nucleotides

    pFF_->ClearGroups();

    if (!isKCAL_)
        conf_data.torsionE *= KJ_TO_KCAL;

    // if not accepted, return
    if (energy_filter[2] < conf_data.torsionE) {
        conf_data.accepted = false;
        return;
    }

    // Get VDW energy
    pFF_->Setup(*current_mol, constraintsTot_); //empty constraint

    conf_data.VDWE = pFF_->E_VDW(false) / n; // Divide by the number of nucleotides

    if (!isKCAL_)
        conf_data.VDWE *= KJ_TO_KCAL;

    // if not accepted, return
    if (energy_filter[3] < conf_data.VDWE || isnan(conf_data.VDWE)) { // nan energies happen when two atoms are in exactly the same position
        conf_data.accepted = false;
        return;
    }

    // Get total energy
    conf_data.total_energy = pFF_->Energy(false) / n;

    if (!isKCAL_)
        conf_data.total_energy *= KJ_TO_KCAL;

    if (energy_filter[4] < conf_data.total_energy || isnan(conf_data.total_energy)) {
        conf_data.accepted = false;
        return;
    }

    // If it passes all energy thresholds, accept backbone candidate
    conf_data.accepted = true;
}

void Chain::setupChain(std::vector<PNAB::Base> &strand, OpenBabel::OBMol &chain, std::vector<unsigned> &new_bond_ids,
                       std::vector<unsigned> &deleted_atoms_ids, std::vector<unsigned> &num_base_unit_atoms,
                       std::vector<unsigned> &bb_start_index, std::vector<double *> &base_coords_vec,
                       std::vector<std::vector<unsigned>> &fixed_bonds_vec, const Backbone &backbone, unsigned chain_index) {

    // Set some variables
    auto bases_A = strand;
    unsigned neighbor_id;
    vector<unsigned> head_ids, tail_ids;
    OBAtom *head, *tail;

    // Get base units
    vector<BaseUnit> bu_A;
    vector<string> base_codes;
    vector<string> base_names;
    for (auto v : bases_A) {
        bu_A.emplace_back(BaseUnit(v, backbone, glycosidic_bond_distance_));
        base_codes.emplace_back(v.getCode());
        base_names.emplace_back(v.getName());
    }

    vector<OBMol> bu_A_mol;
    vector<array<size_t, 2>> bu_linkers_vec;
    vector<vector<vector<unsigned>>> base_fixed_bonds;
    for (auto v : bu_A) {
        // We want a list of linkers for deletion
        bu_linkers_vec.emplace_back(v.getBackboneLinkers());

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

    // Get correct residue number
    int resnum = 0;
    for (int i = 0; i < chain_index; i++) {
        if (build_strand_[i])
            resnum += chain_length_;
    }

    // Create chain by adding all the nucleotides
    auto chain_letter = static_cast<char>('A' + chain_index);
    unsigned c = 0;
    for (auto v : bu_A_mol) {
        for (auto it = v.BeginResidues(); it != v.EndResidues(); ++it) {
            OBResidue *r = *it;
            r->SetName(base_codes[c]);
            r->SetChainNum(chain_index);
            r->SetChain(chain_letter);
            r->SetTitle(base_names[c].c_str());
            if (strand_orientation_[chain_index])
                r->SetNum(resnum + c + 1);
            else
                r->SetNum(resnum + chain_length_ - c);
        }
        chain += v;
        c++;
    }

    // Mark fixed bonds; we want to extract them later
    int num_atoms = 0, n_fixed = 0;
    for (int i=0; i < base_fixed_bonds.size(); i++) {
        for (auto bonds: base_fixed_bonds[i]) {
            OBAtom* a1 = chain.GetAtom(bonds[0] + num_atoms);
            OBAtom* a2 = chain.GetAtom(bonds[1] + num_atoms);
            OBPairData *label1 = new OBPairData;
            OBPairData *label2 = new OBPairData;
            label1->SetAttribute(to_string(n_fixed));
            label2->SetAttribute(to_string(n_fixed));
            a1->CloneData(label1);
            a2->CloneData(label2);
            delete label1;
            delete label2;
            n_fixed++;
        }
        num_atoms += num_base_unit_atoms[i];
    }

    bu_linkers_vec.pop_back();

    std::vector<long int> center_ids;
    unsigned count = 0, index = 0;

    // Determine the connection based on the strand orientation
    if (strand_orientation_[chain_index]) {
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

    // Delete extra hydrogen atoms that remain after connecting the nucleotides
    for (unsigned i = 0; i < head_ids.size(); i++) {
        head = chain.GetAtomById(head_ids[i]);
        tail = chain.GetAtomById(tail_ids[i]);

        // Get indices of the hydrogen atoms
        vector<vector<unsigned long>> hydrogens = {{}, {}};
        FOR_NBORS_OF_ATOM(nbr, tail) {
            if (nbr->GetAtomicNum() == 1)
                hydrogens[0].push_back(nbr->GetId());
        }

        FOR_NBORS_OF_ATOM(nbr, head) {
            if (nbr->GetAtomicNum() == 1)
                hydrogens[1].push_back(nbr->GetId());
        }

        // Delete all the hydrogen atoms from the tail
        for (auto i: hydrogens[0]) {
            deleted_atoms_ids.push_back(i);
            chain.DeleteAtom(chain.GetAtomById(i));
        }

        // Delete hydrogen atoms from the head until you have the correct number of bonds - 1
        // One more bond will be added to the next nucleotide
        for (auto i: hydrogens[1]) {
            if (head->GetExplicitDegree() == head->GetTotalDegree() - 1)
                break;
            deleted_atoms_ids.push_back(i);
            chain.DeleteAtom(chain.GetAtomById(i));
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


    // Get the indices of the fixed bonds after deletion of other atoms
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

    map<unsigned, unsigned> residue_map;
    FOR_RESIDUES_OF_MOL(r, chain) {
        FOR_ATOMS_OF_RESIDUE(a,&*r) {
            residue_map[a->GetIdx()] = r->GetNum();
        }
    }

    vector<int> bond_atoms;

    // Determine the bond and angle atoms only for the first strand
    if (offset == 0) {

        if (chain_length_ != 1) {
            // The bond and angle terms are computed only between the first and second nucleotides
            OBAtom* atom1 = chain.GetAtomById(new_bond_ids[0]);
            OBAtom* atom2 = chain.GetAtomById(new_bond_ids[1]);
            bond_atoms.push_back(atom1->GetIdx());
            bond_atoms.push_back(atom2->GetIdx());

            // The angle atoms are between all triples involving the two linkers and the neighboring atoms
            FOR_NBORS_OF_ATOM(nbr, atom2) {
                if (nbr->GetIdx() != atom1->GetIdx())
                    all_angles_.push_back({atom1->GetIdx() + offset, atom2->GetIdx() + offset, nbr->GetIdx() + offset});
            }
            FOR_NBORS_OF_ATOM(nbr, atom1) {
                if (nbr->GetIdx() != atom2->GetIdx())
                all_angles_.push_back({atom2->GetIdx() + offset, atom1->GetIdx() + offset, nbr->GetIdx() + offset});
            }
        }

        // Ignore all other atoms in the bond calculation
        FOR_ATOMS_OF_MOL(a, chain) {
            if (!(find(bond_atoms.begin(), bond_atoms.end(), a->GetIdx()) != bond_atoms.end()))
                constraintsBond_.AddIgnore(a->GetIdx() + offset);
        }
    }

    // ignore all atoms in the bond calculation for the other strands
    else {
        FOR_ATOMS_OF_MOL(a, chain)
            constraintsBond_.AddIgnore(a->GetIdx() + offset);
    }

    // Determine torsion atoms in the strand
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
                unsigned num = residue_map[nbr->GetIdx()];
                // Exclude dihedral angles in two different residues as these bonds are not
                // rotated in the conformation search
                if ((residue_map[v[1]] == num) &&
                    (residue_map[v[2]] == num) &&
                    (residue_map[nbr2->GetIdx()] == num)) {
                    // These are the correct torsional angles that we want to compute energies for
                    bool fixed = false;
                    // Determine whether this angle is fixed
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

    // Rotation for hexad
    double twist = chain_index * 60.0 * DEG_TO_RAD;
    matrix3x3 z_rot = {{cos(twist), -sin(twist), 0},{sin(twist), cos(twist), 0}, {0, 0, 1}};

    // Reflection for both the duplex and the hexad
    matrix3x3 change_sign = {{1, 0, 0},{0, -1, 0}, {0, 0, -1}};

    //// Reflection that we can apply if we want the same orientation as that of 3DNA
    //matrix3x3 change_sign2 = {{-1, 0, 0},{0, 1, 0}, {0, 0, -1}};

    // Get the global rotation and translations
    auto g_rot = hp.getGlobalRotationMatrix(false, false);
    auto g_trans = hp.getGlobalTranslationVec(false, false);
    unsigned xyzI = 0, local_offset = 0, deleted_atom_index = 0;

    // Get the correct beginning index for the duplex and the hexad
    for (int i=0; i<chain_index; i++) {
        if (build_strand_[i])
            xyzI += 3 * v_chain_[i].NumAtoms();
    }

    // Setup the coordinates for each nucleotide
    for (unsigned i = 0; i < chain_length_; ++i) {
        auto n = num_bu_atoms[i];
        auto r = bb_start_index[i];
        // Get the step translation and rotation 
        auto s_trans = hp.getStepTranslationVec(i, false, false);
        auto s_rot   = hp.getStepRotationMatrix(i, false, false);

        // Set the coordinates for the nucleobases
        double *base_coords = base_coords_vec[i];
        for (unsigned baseI = 0; baseI < r - 1; ++baseI) {
            vector3 v3;
            v3.Set(base_coords + 3 * baseI);

            if (!strand_orientation_[chain_index])
                v3 *= change_sign;

            v3 *= hp.getStepRotationMatrix(1, true, (bool) chain_index);
            v3 *= hp.getGlobalRotationMatrix(true, (bool) chain_index);
            v3 += hp.getStepTranslationVec(1, true, (bool) chain_index);
            v3 += hp.getGlobalTranslationVec(true, (bool) chain_index);
            
            v3 *= g_rot;
            v3 += g_trans;
            v3 *= s_rot;
            v3 += s_trans;

            if (hexad_)
                v3 *= z_rot;

            //// Reflect to get the orientation that agrees with 3DNA for DNA/RNA
            //v3 *=  change_sign2;

            v3.Get(xyz + xyzI);
            xyzI += 3;
        }

        // Set the coordinates for backbone
        unsigned monomer_index = monomer_bb_index_range_[0] - 1;
        for (unsigned bbI = r - 1; bbI < n; bbI++) {
            if (deleted_atoms_ids.empty() || bbI + local_offset + 1 != deleted_atoms_ids[deleted_atom_index] + 1) {
                vector3 v3;
                v3.Set(conf + 3 * monomer_index);

                if (!strand_orientation_[chain_index])
                     v3 *= change_sign;

                v3 *= hp.getStepRotationMatrix(1, true, (bool) chain_index);
                v3 *= hp.getGlobalRotationMatrix(true, (bool) chain_index);
                v3 += hp.getStepTranslationVec(1, true, (bool) chain_index);
                v3 += hp.getGlobalTranslationVec(true, (bool) chain_index);
 
                v3 *= g_rot;
                v3 += g_trans;
                v3 *= s_rot;
                v3 += s_trans;

                if (hexad_)
                    v3 *= z_rot;

                //// Reflect to get the orientation that agrees with 3DNA for DNA/RNA
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
