//
// Created by jbarnett8 on 2/10/18.
//

#include <openbabel/obconversion.h>
#include "Chain.h"

using namespace std;
using namespace PNAB;
using namespace OpenBabel;

Chain::Chain(Bases bases, Backbone backbone, vector<string> strand, string ff_type,
             std::array<unsigned, 2> &range, bool double_stranded, bool hexad) {

    hexad_ = hexad;
    double_stranded_ = double_stranded;
    ff_type_ = ff_type;
    monomer_bb_index_range_ = range;
    chain_length_ = static_cast<unsigned>(strand.size());
    backbone_ = backbone;
    chain_length_ = static_cast<unsigned>(strand.size());
    pFF_ = OBForceField::FindForceField(ff_type_);
    isKCAL_ = pFF_->GetUnit().find("kcal") != string::npos;

    if (double_stranded_) {
        auto bases_A = bases.getBasesFromStrand(strand);
        auto bases_B = bases.getComplimentBasesFromStrand(strand);
        if (bases_B.empty()) {
            cerr << "User requested double stranded chain, however not all bases have a valid compliment. Please check"
                 << " input file. Exiting..." << endl;
            throw 1;
        }

        // Strand 1
        setupChain(bases_A, chain_, new_bond_ids_, deleted_atoms_ids_, num_bu_A_mol_atoms_, bb_start_index_,
                   base_coords_vec_, 0, double_stranded_, hexad);
        std::sort(deleted_atoms_ids_.begin(), deleted_atoms_ids_.end());
        setupFFConstraints(chain_, new_bond_ids_);
        combined_chain_ += chain_;

        // Strand 2
        setupChain(bases_B, c_chain_, c_new_bond_ids_, c_deleted_atoms_id_, num_bu_B_mol_atoms_, c_bb_start_index_,
                   c_base_coords_vec_, 1, double_stranded_, false);
        std::sort(c_deleted_atoms_id_.begin(), c_deleted_atoms_id_.end());
        setupFFConstraints(c_chain_, c_new_bond_ids_, c_chain_.NumAtoms());
        combined_chain_ += c_chain_;
    }

    else if (hexad_) {
        auto bases_A = bases.getBasesFromStrand(strand);
        auto bases_B = bases.getComplimentBasesFromStrand(strand);
        if (bases_B.empty()) {
            cerr << "User requested hexad chain, however not all bases have a valid compliment. Please check"
                 << " input file. Exiting..." << endl;
            throw 1;
        }

        // Strand 1
        setupChain(bases_A, chain_, new_bond_ids_, deleted_atoms_ids_, num_bu_A_mol_atoms_, bb_start_index_,
                   base_coords_vec_, 0, false, hexad);
        std::sort(deleted_atoms_ids_.begin(), deleted_atoms_ids_.end());
        setupFFConstraints(chain_, new_bond_ids_);
        combined_chain_ += chain_;

        // Strand 2
        setupChain(bases_B, c_chain_, c_new_bond_ids_, c_deleted_atoms_id_, num_bu_B_mol_atoms_, c_bb_start_index_,
                   c_base_coords_vec_, 1, false, hexad_);
        std::sort(c_deleted_atoms_id_.begin(), c_deleted_atoms_id_.end());
        setupFFConstraints(c_chain_, c_new_bond_ids_, chain_.NumAtoms());
        combined_chain_ += c_chain_;


        // Strand 3
        setupChain(bases_A, d_chain_, d_new_bond_ids_, d_deleted_atoms_id_, num_bu_A_mol_atoms_, d_bb_start_index_,
                   d_base_coords_vec_, 2, false, hexad_);
        std::sort(d_deleted_atoms_id_.begin(), d_deleted_atoms_id_.end());
        setupFFConstraints(d_chain_, d_new_bond_ids_, chain_.NumAtoms() + c_chain_.NumAtoms());
        combined_chain_ += d_chain_;

        // Strand 4
        setupChain(bases_B, e_chain_, e_new_bond_ids_, e_deleted_atoms_id_, num_bu_B_mol_atoms_, e_bb_start_index_,
                   e_base_coords_vec_, 3, false, hexad_);
        std::sort(e_deleted_atoms_id_.begin(), e_deleted_atoms_id_.end());
        setupFFConstraints(e_chain_, e_new_bond_ids_, chain_.NumAtoms() * 2 + c_chain_.NumAtoms());
        combined_chain_ += e_chain_;

        // Strand 5
        setupChain(bases_A, f_chain_, f_new_bond_ids_, f_deleted_atoms_id_, num_bu_A_mol_atoms_, f_bb_start_index_,
                   f_base_coords_vec_, 4, false, hexad_);
        std::sort(f_deleted_atoms_id_.begin(), f_deleted_atoms_id_.end());
        setupFFConstraints(f_chain_, f_new_bond_ids_, chain_.NumAtoms() * 2 + c_chain_.NumAtoms() * 2);
        combined_chain_ += f_chain_;

        // Strand 6
        setupChain(bases_B, g_chain_, g_new_bond_ids_, g_deleted_atoms_id_, num_bu_B_mol_atoms_, g_bb_start_index_,
                   g_base_coords_vec_, 5, false, hexad_);
        std::sort(g_deleted_atoms_id_.begin(), g_deleted_atoms_id_.end());
        setupFFConstraints(g_chain_, g_new_bond_ids_, chain_.NumAtoms() * 3 + c_chain_.NumAtoms() * 2);
        combined_chain_ += g_chain_;

    }    


    else {
        auto bases_A = bases.getBasesFromStrand(strand);
        setupChain(bases_A, chain_, new_bond_ids_, deleted_atoms_ids_, num_bu_A_mol_atoms_, bb_start_index_,
                   base_coords_vec_, 0, double_stranded_, hexad);
        std::sort(deleted_atoms_ids_.begin(), deleted_atoms_ids_.end());
        setupFFConstraints(chain_, new_bond_ids_);
    }

    if (double_stranded_ || hexad_)
        pFF_->Setup(combined_chain_);
    else
        pFF_->Setup(chain_);
}

ConformerData Chain::generateConformerData(double *conf, HelicalParameters &hp) {
    if (chain_.NumAtoms() <= 0) {
        cout << "There are no atoms in the chain_. Exiting..." << endl;
        exit(0);
    }

    unsigned num_cooords;
    if (double_stranded_)
        num_cooords = (chain_.NumAtoms() + c_chain_.NumAtoms()) * 3;
    else if (hexad_)
        num_cooords = (chain_.NumAtoms() * 3 + c_chain_.NumAtoms() * 3) * 3;
    else
        num_cooords = chain_.NumAtoms() * 3;
    auto *xyz = new double[num_cooords];


    if (double_stranded_) {
        setCoordsForChain(xyz,conf,hp,num_bu_A_mol_atoms_,bb_start_index_,base_coords_vec_,deleted_atoms_ids_, false, false, 0);
        setCoordsForChain(xyz,conf,hp,num_bu_B_mol_atoms_,c_bb_start_index_,c_base_coords_vec_,c_deleted_atoms_id_, true, false, 1);
    }

    else if (hexad_) {
        setCoordsForChain(xyz,conf,hp,num_bu_A_mol_atoms_,bb_start_index_,base_coords_vec_,deleted_atoms_ids_, false, false, 0);
        setCoordsForChain(xyz,conf,hp,num_bu_B_mol_atoms_,c_bb_start_index_,c_base_coords_vec_,c_deleted_atoms_id_, false, true, 1);
        setCoordsForChain(xyz,conf,hp,num_bu_A_mol_atoms_,d_bb_start_index_,d_base_coords_vec_,d_deleted_atoms_id_, false, true, 2);
        setCoordsForChain(xyz,conf,hp,num_bu_B_mol_atoms_,e_bb_start_index_,e_base_coords_vec_,e_deleted_atoms_id_, false, true, 3);
        setCoordsForChain(xyz,conf,hp,num_bu_A_mol_atoms_,f_bb_start_index_,f_base_coords_vec_,f_deleted_atoms_id_, false, true, 4);
        setCoordsForChain(xyz,conf,hp,num_bu_B_mol_atoms_,g_bb_start_index_,g_base_coords_vec_,g_deleted_atoms_id_, false, true, 5);
    }
    else {
        setCoordsForChain(xyz,conf,hp,num_bu_A_mol_atoms_,bb_start_index_,base_coords_vec_,deleted_atoms_ids_, false, false, 0);
    }

    ConformerData data;
    data.coords = xyz;
    data.chain_coords_present = true;
    fillConformerEnergyData(xyz, data);
    return data;
}

void Chain::fillConformerEnergyData(double *xyz, PNAB::ConformerData &conf_data) {

    OBMol *current_mol;

    if (double_stranded_ || hexad_)
        current_mol = &combined_chain_;
    else
        current_mol = &chain_;

    current_mol->SetCoordinates(xyz);

    auto n = chain_length_;
    if (double_stranded_)
        n *= 2;
    else if (hexad_)
        n *= 6;

    // Get total energies and VDW
    pFF_->Setup(*current_mol, constraintsTot_);
    conf_data.total_energy = pFF_->Energy() / n;
    conf_data.VDWE = pFF_->E_VDW() / n;
    conf_data.totTorsionE = pFF_->E_Torsion() / n;

    if (n - 1 <= 0)
        n = 2;

    // Get torsional energies
    pFF_->SetConstraints(constraintsTor_);
    conf_data.torsionE = pFF_->E_Torsion() / (n - 1);

    // Get angle energy
    pFF_->SetConstraints(constraintsAng_);
    conf_data.angleE = pFF_->E_Angle() / (n - 1);

    // Get bond energy
    pFF_->SetConstraints(constraintsBond_);
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

void Chain::setupChain(std::vector<PNAB::Base> &strand, OpenBabel::OBMol &chain, std::vector<unsigned> &new_bond_ids,
                       std::vector<unsigned> &deleted_atoms_ids, std::vector<unsigned> &num_base_unit_atoms,
                       std::vector<unsigned> &bb_start_index, std::vector<double *> &base_coords_vec,
                       unsigned chain_index, bool double_stranded, bool hexad) {

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
            if (chain_index % 2 == 0 || hexad) {
                r->SetNum(c + 1);
            }
            else if (chain_index % 2 == 1) {
                r->SetNum(bu_A_mol.size()*2 -  c);
            }
        }
        chain += v;
        c++;
    }

    bu_linkers_vec.pop_back();

    std::vector<long int> center_ids;
    unsigned count = 0, index = 0;
    if (chain_index == 0 || hexad) {
        for (unsigned i = 0; i < bu_linkers_vec.size(); ++i) {
            head_ids.push_back(static_cast<int>(chain.GetAtom(count + bu_linkers_vec[i][0])->GetId()));
            count += num_base_unit_atoms[index++];
            tail_ids.push_back(static_cast<int>(chain.GetAtom(count + bu_linkers_vec[i + 1][1])->GetId()));
        }
    }
    else if (chain_index == 1)
        for (unsigned i = 0; i < bu_linkers_vec.size(); ++i) {
            tail_ids.push_back(static_cast<int>(chain.GetAtom(count + bu_linkers_vec[i    ][1])->GetId()));
            count += num_base_unit_atoms[index++];
            head_ids.push_back(static_cast<int>(chain.GetAtom(count + bu_linkers_vec[i + 1][0])->GetId()));
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
            constraintsTor_.AddIgnore(a->GetIdx() + offset);
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

}

void Chain::setCoordsForChain(double *xyz, double *conf, PNAB::HelicalParameters &hp,
                              std::vector<unsigned> &num_bu_atoms,
                              std::vector<unsigned> &bb_start_index, std::vector<double *> &base_coords_vec,
                              std::vector<unsigned> &deleted_atoms_ids, bool is_second_strand, bool is_hexad, unsigned chain_index) {

    vector3 z_trans{0, 0, 0};
    vector3 s_trans_prev, g_trans_prev;
    matrix3x3 s_rot_prev, g_rot_prev, z_rot;
    auto g_rot = hp.getGlobalRotationOBMatrix();
    auto g_trans = hp.getGlobalTranslationVec();
    unsigned xyzI = 0, local_offset = 0, deleted_atom_index = 0;
    if (is_second_strand)
        xyzI = 3 * chain_.NumAtoms();
    else if (is_hexad)
        xyzI = 3 * chain_.NumAtoms() * ((chain_index + 1) / 2) + 3 * c_chain_.NumAtoms() * (chain_index / 2); 


    matrix3x3 change_sign({1, 0, 0},{0, -1, 0}, {0, 0, -1});
    if (is_hexad) {
        double twist = chain_index * 60.0 * DEG_TO_RAD;
        z_rot = matrix3x3{{cos(twist), -sin(twist), 0},{sin(twist), cos(twist), 0}, {0, 0, 1}};
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

            if (is_second_strand)
                v3 *= change_sign;

            else if (is_hexad) {
                v3 *= z_rot;
            }

            v3 += s_trans; v3 *= s_rot;
            v3.Get(xyz + xyzI);
            xyzI += 3;
        }
        unsigned monomer_index = monomer_bb_index_range_[0] - 1;
        for (unsigned bbI = r - 1; bbI < n; bbI++) {
            if (deleted_atoms_ids.empty() || bbI + local_offset + 1 != deleted_atoms_ids[deleted_atom_index] + 1) {
                vector3 v3;
                v3.Set(conf + 3 * monomer_index);

                v3 += g_trans; v3 *= g_rot;

                if (is_second_strand)
                    v3 *= change_sign;

                else if (is_hexad) {
                    v3 *= z_rot;
                }
                
                v3 += s_trans; v3 *= s_rot;
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
