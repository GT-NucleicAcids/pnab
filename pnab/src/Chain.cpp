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
                   v_base_coords_vec_[i], i);
        std::sort(v_deleted_atoms_ids_[i].begin(), v_deleted_atoms_ids_[i].end());
        setupFFConstraints(v_chain_[i], v_new_bond_ids_[i], offset);
        offset += v_chain_[i].NumAtoms();
        combined_chain_ += v_chain_[i];
    }

    pFF_->Setup(combined_chain_);
    checkTorsions();
}


void Chain::checkTorsions() {
    // Get whether we should read the torsion
    pFF_->GetAtomTypes(combined_chain_);
    ostringstream torsion_out;
    pFF_->SetLogFile(&torsion_out);
    pFF_->SetLogLevel(3);
    pFF_->SetConstraints(constraintsTor_);
    double torsionE = pFF_->E_Torsion(false);
    pFF_->SetLogLevel(0);

    string torsion_out_str = torsion_out.str();
    std::istringstream string_stream(torsion_out_str);
    string line;
    int i = -1, j=0, angle_position;
    while(getline(string_stream, line, '\n')) {
        if (j==0) {
            size_t found = line.find("ATOM TYPES");
            if (found!=std::string::npos) {
                j++;
                continue;
            }
        }

        stringstream line_string_stream(line);
        string token;
        vector<string> components;
        while (line_string_stream >> token) {
            components.push_back(token);
        }

        if (j==1) {
            for (int k=0; k < components.size(); k++) {
                if (components[k].find("ANGLE") !=std::string::npos)
                    angle_position = k;
            }
            j++;
            continue;
        }
        if (j==2) {
            i++;
            if (i == number_of_all_torsions)
                break;
            bool save = true;
            for (int k = 0; k < all_torsions_.size(); k++) {
                vector<unsigned int> v = all_torsions_[k];
                auto a1 = combined_chain_.GetAtom((unsigned int) all_torsions_[k][0]);
                auto a2 = combined_chain_.GetAtom((unsigned int) all_torsions_[k][1]);
                auto a3 = combined_chain_.GetAtom((unsigned int) all_torsions_[k][2]);
                auto a4 = combined_chain_.GetAtom((unsigned int) all_torsions_[k][3]);
                OBPairData *t1 = (OBPairData*) a1->GetData("FFAtomType");
                OBPairData *t2 = (OBPairData*) a2->GetData("FFAtomType");
                OBPairData *t3 = (OBPairData*) a3->GetData("FFAtomType");
                OBPairData *t4 = (OBPairData*) a4->GetData("FFAtomType");
                int t = (int) (1000.0 * combined_chain_.GetTorsion(a1, a2, a3, a4));

                if (((components[0] == t1->GetValue() && components[1] == t2->GetValue() && components[2] == t3->GetValue() && components[3] == t4->GetValue()) || 
                   (components[3] == t1->GetValue() && components[2] == t2->GetValue() && components[1] == t3->GetValue() && components[0] == t4->GetValue()))
                    && abs(t - (int) (1000.0 * atof(components[angle_position].c_str())))) {
                    correct_torsions.push_back(true);
                    save = false;
                    break;
                }

            }

            if (save) {
                correct_torsions.push_back(false);
            }

        }
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

    pFF_->Setup(*current_mol);

    if (!(n == 1) || !(double_stranded_ && n == 2) || !(hexad_ && n==6)) {

        // Get bond energy
        pFF_->SetConstraints(constraintsBond_);
        conf_data.bondE = pFF_->E_Bond(false) / (n - subtract);
        if (!isKCAL_)
            conf_data.bondE *= KJ_TO_KCAL;

        if (energy_filter[0] < conf_data.bondE) {
            conf_data.accepted = false;
            return;
        }


        // Get angle energy
        pFF_->SetConstraints(constraintsAng_);
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
    // Read correct torsions from Openbabel logging
    ostringstream torsion_out;
    pFF_->SetLogFile(&torsion_out);
    pFF_->SetLogLevel(3);
    pFF_->SetConstraints(constraintsTor_);
    pFF_->E_Torsion(false);
    pFF_->SetLogLevel(0);

    double E_tor = 0.0;
    string torsion_out_str = torsion_out.str();
    std::istringstream string_stream(torsion_out_str);
    string line;
    bool read = false;
    int i = -1;
    while(getline(string_stream, line, '\n')) {
        size_t found = line.find("---------------------");
        if (found!=std::string::npos) {
            read = true;
            continue;
        }

        stringstream line_string_stream(line);
        string token;
        vector<string> components;
        if (read) {
            i++;
            if (i == number_of_all_torsions)
                break;
            if (!correct_torsions[i])
                continue;
            while (line_string_stream >> token) {
                components.push_back(token);
            }
            E_tor += atof(components[components.size()-1].c_str());
        }
    }

    conf_data.torsionE = E_tor / n ;

    if (!isKCAL_)
        conf_data.torsionE *= KJ_TO_KCAL;

    if (energy_filter[2] < conf_data.torsionE) {
        conf_data.accepted = false;
        return;
    }


    OBFFConstraints constraints;
    pFF_->SetConstraints(constraints);

    // Get VDW energy
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
                       unsigned chain_index) {

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

};

void Chain::setupFFConstraints(OpenBabel::OBMol &chain, std::vector<unsigned> &new_bond_ids, unsigned offset) {
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

    set<int> atoms;

    while(r) {
        vector<int> v = r->GetDihedralAtoms();
        FOR_NBORS_OF_ATOM(nbr, chain.GetAtom(v[1])) {
            if (nbr->GetIdx() == v[2])
                continue;
            FOR_NBORS_OF_ATOM(nbr2, chain.GetAtom(v[2])) {
                if (nbr2->GetIdx() == v[1])
                    continue;
                // These are the correct torsional angles that we want to compute energies for
                all_torsions_.push_back({nbr->GetIdx() + offset, v[1] + offset, v[2] + offset, nbr2->GetIdx() + offset});
                atoms.insert(nbr->GetIdx());
                atoms.insert(v[1]);
                atoms.insert(v[2]);
                atoms.insert(nbr2->GetIdx());
            }
        }

        r = rl.NextRotor(ri);
    }

    // Ignore all atoms not involved in the torsion
    // This might count additional torsional angles that are incorrect. This happens for
    // angles that involve nonrotatable bonds in the molecule but are not ignored, i.e. terminal atoms in rotatable bonds.
    FOR_ATOMS_OF_MOL(a, chain) {
        if (!(atoms.find(a->GetIdx()) != atoms.end()))
            constraintsTor_.AddIgnore(a->GetIdx() + offset);
    }

    // Here we determine the number of all torsional angles that openbabel would calculate energies for, correct or otherwise.
    FOR_TORSIONS_OF_MOL(t, chain) {
        int a = chain.GetAtom((*t)[0] + 1)->GetIdx();
        int b = chain.GetAtom((*t)[1] + 1)->GetIdx();
        int c = chain.GetAtom((*t)[2] + 1)->GetIdx();
        int d = chain.GetAtom((*t)[3] + 1)->GetIdx();
        if ((atoms.find(a) != atoms.end()) && (atoms.find(b) != atoms.end()) && (atoms.find(c) != atoms.end()) && (atoms.find(d) != atoms.end())) {
            bool save = true;
            for (auto v: all_torsions_) {
                if (((a == v[0] && b == v[1] && c == v[2] && d == v[3]) || (a == v[3] && b == v[2] && c == v[1] && d == v[0]))) {
                    save = false;
                    number_of_all_torsions++;
                    break;
                }
            }
            if (save) {
                number_of_all_torsions++;
            }
        }
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
