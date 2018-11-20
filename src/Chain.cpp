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

    double_stranded_ = double_stranded;
    ff_type_ = ff_type;
    monomer_bb_index_range_ = range;
    chain_length_ = static_cast<unsigned>(strand.size());
    backbone_ = backbone;
    chain_length_ = static_cast<unsigned>(strand.size());
    pFF_ = OBForceField::FindForceField(ff_type_);
    isKCAL_ = pFF_->GetUnit().find("kcal") != string::npos;

    // Setting up constraints
    constraintsTot_ = new OBFFConstraints;
    constraintsTor_ = new OBFFConstraints;
    constraintsAng_ = new OBFFConstraints;
    constraintsBond_ = new OBFFConstraints;

    auto bases_A = bases.getBasesFromStrand(strand);

    setupChain(bases_A, chain_, new_bond_ids_, deleted_atoms_ids_, num_bu_A_mol_atoms_, bb_start_index_,
               base_coords_vec_, n1_or_n3_tuple_, 0);
    std::sort(deleted_atoms_ids_.begin(), deleted_atoms_ids_.end());

    setupFFConstraints(chain_, new_bond_ids_);

    if (double_stranded_) {
        auto bases_B = bases.getComplimentBasesFromStrand(strand);
        if (bases_B.empty()) {
            cerr << "User requested double stranded chain, however not all bases have a valid compliment. Please check"
                 << " input file. Exiting..." << endl;
            throw 1;
        }
        reverse(bases_B.begin(),bases_B.end());
        setupChain(bases_B, c_chain_, c_new_bond_ids_, c_deleted_atoms_id_, num_bu_B_mol_atoms_, c_bb_start_index_,
                   c_base_coords_vec_, c_n1_or_n3_tuple_, 1);
        std::sort(c_deleted_atoms_id_.begin(), c_deleted_atoms_id_.end());

        setupFFConstraints(c_chain_, c_new_bond_ids_, chain_.NumAtoms());

        combined_chain_ += chain_;
        combined_chain_ += c_chain_;
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

    if (double_stranded_)
        pFF_->Setup(combined_chain_);
    else
        pFF_->Setup(chain_);
}

ConformerData Chain::generateConformerData(double *conf, HelicalParameters &hp) {
    if (chain_.NumAtoms() <= 0) {
        cout << "There are no atoms in the chain_. Exiting..." << endl;
        exit(0);
    }

    // new new *new* code
    auto num_cooords = double_stranded_ ? (chain_.NumAtoms() + c_chain_.NumAtoms()) * 3 : chain_.NumAtoms() * 3;
    auto *xyz = new double[num_cooords];

    setCoordsForChain(xyz,conf,hp,num_bu_A_mol_atoms_,bb_start_index_,base_coords_vec_,deleted_atoms_ids_, false);

    if (double_stranded_) {
        setCoordsForChain(xyz,conf,hp,num_bu_B_mol_atoms_,c_bb_start_index_,c_base_coords_vec_,c_deleted_atoms_id_, true);
    }

    ConformerData data;
    data.coords = xyz;
    data.chain_coords_present = true;
    fillConformerEnergyData(xyz, data);
    return data;
}

void Chain::fillConformerEnergyData(double *xyz, PNAB::ConformerData &conf_data) {

    OBMol *current_mol;

    if (double_stranded_)
        current_mol = &combined_chain_;
    else
        current_mol = &chain_;

    current_mol->SetCoordinates(xyz);

    auto n = chain_length_;

    // Get total energies and VDW
    pFF_->Setup(*current_mol, *constraintsTot_);
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

void Chain::setupChain(std::vector<PNAB::Base> &strand, OpenBabel::OBMol &chain, std::vector<unsigned> &new_bond_ids,
                       std::vector<unsigned> &deleted_atoms_ids, std::vector<unsigned> &num_base_unit_atoms,
                       std::vector<unsigned> &bb_start_index, std::vector<double *> &base_coords_vec,
                       std::vector<std::tuple<std::string, unsigned, unsigned>> &h1_or_h3_tuple, unsigned chain_index) {
    if (chain_index > 1) {
        cerr << "Program only supports up to two chains. This error really shouldn't trigger." << endl;
        throw 1;
    }

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

        if (double_stranded_) {
            h1_or_h3_tuple.emplace_back(v.getC6AndN3AtomIndices());
        }

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
    //TODO fix this inconsistency temporarily added for quick fix
    if (chain_index == 0 || chain_index == 1) {
        if (double_stranded_) {
            unsigned n_count = 0;
            for (unsigned i = 0; i < chain_length_; ++i) {
                get<1>(h1_or_h3_tuple[i]) = chain.GetAtom(n_count + get<1>(h1_or_h3_tuple[i]))->GetId();
                if (get<0>(h1_or_h3_tuple[i]).find("N3") != string::npos)
                    get<2>(h1_or_h3_tuple[i]) = chain.GetAtom(n_count + get<2>(h1_or_h3_tuple[i]))->GetId();
                n_count += num_base_unit_atoms[i];
            }
        }
        for (unsigned i = 0; i < bu_linkers_vec.size(); ++i) {
            head_ids.push_back(static_cast<int>(chain.GetAtom(count + bu_linkers_vec[i][0])->GetId()));
            count += num_base_unit_atoms[index++];
            tail_ids.push_back(static_cast<int>(chain.GetAtom(count + bu_linkers_vec[i + 1][1])->GetId()));
        }
    }  else if (chain_index == 1)
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
        if (chain_index == 0) {
            FOR_NBORS_OF_ATOM(nbr, tail) {
                neighbor_id = nbr->GetId();
                break;
            }
            deleted_atoms_ids.push_back(tail->GetId());
            chain.DeleteAtom(tail);
            new_bond_ids.push_back(head->GetId());
            new_bond_ids.push_back(neighbor_id);
            chain.AddBond(head->GetIdx(), chain.GetAtomById(neighbor_id)->GetIdx(), 1);
        } else if (chain_index == 1) {
            FOR_NBORS_OF_ATOM(nbr, head) {
                neighbor_id = nbr->GetId();
                break;
            }
            deleted_atoms_ids.push_back(head->GetId());
            chain.DeleteAtom(head);
            new_bond_ids.push_back(tail->GetId());
            new_bond_ids.push_back(neighbor_id);
            chain.AddBond(tail->GetIdx(), chain.GetAtomById(neighbor_id)->GetIdx(), 1);
        }
    }

    if (double_stranded_) {
        for (auto &v : h1_or_h3_tuple) {
            get<1>(v) = chain.GetAtomById(get<1>(v))->GetIdx();
            if (get<0>(v).find("N3") != string::npos)
                get<2>(v) = chain.GetAtomById(get<2>(v))->GetIdx();
        }


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
            constraintsTor_->AddIgnore(a->GetIdx() + offset);
            constraintsAng_->AddIgnore(a->GetIdx() + offset);
            constraintsBond_->AddIgnore(a->GetIdx() + offset);
        }
    }

    // Angle energy, must fix non-angle atoms
    for (int i = 1; i < new_bond_ids.size(); i+=2) {
        FOR_NBORS_OF_ATOM(nbr, chain.GetAtomById(new_bond_ids.at(i))) {
            nbrIdx = nbr->GetIdx();
            if (nbrIdx != chain.GetAtomById(new_bond_ids.at(i - 1))->GetIdx()) {
                constraintsAng_->AddIgnore(nbrIdx + offset);
                constraintsBond_->AddIgnore(nbrIdx + offset);
            }
        }
    }

    // Bond energy, must fix non-bond atoms
    for (int i = 0; i < new_bond_ids.size(); i+=2) {
        FOR_NBORS_OF_ATOM(nbr, chain.GetAtomById(new_bond_ids.at(i))) {
            nbrIdx = nbr->GetIdx();
            if (nbrIdx != chain.GetAtomById(new_bond_ids.at(i + 1))->GetIdx()) {
                constraintsBond_->AddIgnore(nbrIdx + offset);
            }
        }
    }

    if (!(constraintsTot_ && constraintsTor_ && constraintsAng_ && constraintsBond_)) {
        cerr << "Error initializing force field constraints for chain. Exiting..." << endl;
        exit(1);
    }
}

void Chain::setCoordsForChain(double *xyz, double *conf, PNAB::HelicalParameters &hp,
                              std::vector<unsigned> &num_bu_atoms,
                              std::vector<unsigned> &bb_start_index, std::vector<double *> &base_coords_vec,
                              std::vector<unsigned> &deleted_atoms_ids, bool is_second_strand) {

    vector3 z_trans{0, 0, 0};
    vector3 s_trans_prev, g_trans_prev;
    matrix3x3 s_rot_prev, g_rot_prev;
    auto g_rot = hp.getGlobalRotationOBMatrix();
    auto g_trans = hp.getGlobalTranslationVec();
    unsigned xyzI = is_second_strand ? 3 * chain_.NumAtoms() : 0,
             local_offset = 0, deleted_atom_index = 0;

    for (unsigned i = 0; i < chain_length_; ++i) {
        auto n = num_bu_atoms[i];
        auto r = bb_start_index[i];
        auto s_trans = hp.getStepTranslationVec(i);
        auto s_rot   = hp.getStepRotationOBMatrix(i);
        double *base_coords = base_coords_vec[i];
        for (unsigned baseI = 0; baseI < r - 1; ++baseI) {
            vector3 v3;
            v3.Set(base_coords + 3 * baseI);
//            if (is_second_strand)
//                v3 *= change_sign;
            v3 += g_trans; v3 *= g_rot; v3 += s_trans; v3 *= s_rot;
            v3.Get(xyz + xyzI);
            xyzI += 3;
        }
        unsigned monomer_index = monomer_bb_index_range_[0] - 1;
        for (unsigned bbI = r - 1; bbI < n; bbI++) {
            if (deleted_atoms_ids.empty() || bbI + local_offset + 1 != deleted_atoms_ids[deleted_atom_index] + 1) {
                vector3 v3;
                v3.Set(conf + 3 * monomer_index);
//                if (is_second_strand)
//                    v3 *= change_sign;
                v3 += g_trans; v3 *= g_rot; v3 += s_trans; v3 *= s_rot;
                v3.Get(xyz + xyzI);
                xyzI += 3;
            } else {
                deleted_atom_index++;
            }
            monomer_index++;
        }
        local_offset += n;
    }
//    if (double_stranded_ && !is_second_strand) {
//        auto l = chain_.NumAtoms();
//
//        center_ = {0, 0, 0};
//        for (unsigned i = 0; i < l; ++i) {
//            auto z = mol_weights_[i] * xyz[3 * i + 2];
////            auto z = xyz[3 * i + 2] / l;
////            auto y = xyz[3 * i + 1] / l;
//            //auto x = xyz[3 * i    ] / l;
//            center_ += {0, 0, z};
//        }
////        center_ = {0, 0, xyz[(center_indices_[0] - 1) * 3 + 2]};
////        if (center_indices_.size() == 2) {
////            center_ += {0, 0, xyz[(center_indices_[1] - 1) * 3 + 2]};
////            center_ = center_ / 2;
////        }
//        for (unsigned i = 0; i < l; ++i) {
//            vector3 v3;
//            v3.Set(xyz + 3 * i);
//            v3 -= center_;
//            v3.Get(xyz + 3 * i);
//        }
//    }
    if (double_stranded_ && is_second_strand) {

        auto shift_index = 3 * chain_.NumAtoms();
        auto l = c_chain_.NumAtoms();


        //        double twist = 36.5 * DEG_TO_RAD; //dna
        double twist = 33.2 * DEG_TO_RAD;   //rna


//        double twist = hp.getTwist() * DEG_TO_RAD;
        unsigned n = chain_length_ - 1;
        auto z_rot = matrix3x3{{cos(twist * n), -sin(twist * n), 0},{sin(twist * n), cos(twist * n), 0}, {0, 0, 1}};








        if (!z_trans.CanBeNormalized()) {
            z_trans = {0, 0, 0};
            std::reverse(c_n1_or_n3_tuple_.begin(), c_n1_or_n3_tuple_.end());
            unsigned long range = n1_or_n3_tuple_.size();
            for (unsigned tuple_i = 0; tuple_i < range; ++tuple_i) {
                auto first_base = n1_or_n3_tuple_[tuple_i];
                auto last_c_base = c_n1_or_n3_tuple_[tuple_i];
                auto i = get<1>(first_base) - 1, j = get<1>(last_c_base) - 1;
                vector3 c6_vec, n1_vec, n3_vec, n1_n3_vec;
                if (get<0>(first_base).find("N3") != string::npos) {
                    // i refers to N3 and j refers to N1
                    // cout << "first base" << endl;
                    auto c6_index = get<2>(first_base) - 1;
                    c6_vec = vector3{xyz[3 * c6_index], xyz[3 * c6_index + 1], xyz[3 * c6_index + 2]};
                    n3_vec = vector3(xyz[3 * i], xyz[3 * i + 1], xyz[3 * i + 2]);
                    n1_vec = vector3(xyz[3 * j + shift_index], xyz[3 * j + 1 + shift_index],
                                     xyz[3 * j + 2 + shift_index]);

                    if (tuple_i == 0)
                        twist = calcRotation(0, 1E-12, 3600, c6_vec, n1_vec, n3_vec);
                    c6_vec = n3_vec - c6_vec;
                    n1_n3_vec = n1_vec - n3_vec;

                } else if (get<0>(last_c_base).find("N3") != string::npos) {

                    // i refers to N1 and j refers to N3
                    // cout << "lase base" << endl;
                    auto c6_index = get<2>(last_c_base) - 1;
                    c6_vec = vector3{xyz[3 * c6_index + shift_index], xyz[3 * c6_index + 1 + shift_index],
                                     xyz[3 * c6_index + 2 + shift_index]};
                    n1_vec = vector3(xyz[3 * i], xyz[3 * i + 1], xyz[3 * i + 2]);
                    n3_vec = vector3(xyz[3 * j + shift_index], xyz[3 * j + 1 + shift_index],
                                     xyz[3 * j + 2 + shift_index]);
                    c6_vec = n3_vec - c6_vec;
                    n1_n3_vec = n3_vec - n1_vec;

                } else {
                    cerr << "There was an error making the second strand. Compliment bases must match an N3 to N1."
                         << endl;
                    throw "NoN3ToN1Exception";
                }
            }
            // TODO make this more efficient
            // Reverse the tuples back, probably not necessary
            std::reverse(c_n1_or_n3_tuple_.begin(), c_n1_or_n3_tuple_.end());
        }


        z_rot = matrix3x3{{cos(twist * n), -sin(twist * n), 0},{sin(twist * n), cos(twist * n), 0}, {0, 0, 1}};






        matrix3x3 change_sign({1, 0, 0},{0, -1, 0}, {0, 0, -1});
        for (unsigned i = 0; i < l; ++i) {
            vector3 v3;
            v3.Set(xyz + 3 * i + shift_index);
            v3 *= change_sign;
            v3 *= z_rot;
            v3.Get(xyz + 3 * i + shift_index);
        }




        if (!z_trans.CanBeNormalized()) {
            z_trans = {0, 0, 0};
            std::reverse(c_n1_or_n3_tuple_.begin(), c_n1_or_n3_tuple_.end());
            unsigned long range = n1_or_n3_tuple_.size();
            for (unsigned tuple_i = 0; tuple_i < range; ++tuple_i) {
                auto first_base = n1_or_n3_tuple_[tuple_i];
                auto last_c_base = c_n1_or_n3_tuple_[tuple_i];
                auto i = get<1>(first_base) - 1, j = get<1>(last_c_base) - 1;
                vector3 c6_vec, n1_vec, n3_vec, n1_n3_vec;
                if (get<0>(first_base).find("N3") != string::npos) {
                    // i refers to N3 and j refers to N1
                    // cout << "first base" << endl;
                    auto c6_index = get<2>(first_base) - 1;
                    c6_vec = vector3{xyz[3 * c6_index], xyz[3 * c6_index + 1], xyz[3 * c6_index + 2]};
                    n3_vec = vector3(xyz[3 * i], xyz[3 * i + 1], xyz[3 * i + 2]);
                    n1_vec = vector3(xyz[3 * j + shift_index], xyz[3 * j + 1 + shift_index],
                                     xyz[3 * j + 2 + shift_index]);

//                    if (tuple_i == 0)
//                        twist = calcRotation(0, 1E-12, 3600, c6_vec, n1_vec, n3_vec);
                    c6_vec = n3_vec - c6_vec;
                    n1_n3_vec = n1_vec - n3_vec;

                } else if (get<0>(last_c_base).find("N3") != string::npos) {

                    // i refers to N1 and j refers to N3
                    // cout << "lase base" << endl;
                    auto c6_index = get<2>(last_c_base) - 1;
                    c6_vec = vector3{xyz[3 * c6_index + shift_index], xyz[3 * c6_index + 1 + shift_index],
                                          xyz[3 * c6_index + 2 + shift_index]};
                    n1_vec = vector3(xyz[3 * i], xyz[3 * i + 1], xyz[3 * i + 2]);
                    n3_vec = vector3(xyz[3 * j + shift_index], xyz[3 * j + 1 + shift_index],
                                     xyz[3 * j + 2 + shift_index]);
                    c6_vec = n3_vec - c6_vec;
                    n1_n3_vec = n3_vec - n1_vec;
//                                     vector3(xyz[3 * i], xyz[3 * i + 1], xyz[3 * i + 2]);

                } else {
                    cerr << "There was an error making the second strand. Compliment bases must match an N3 to N1."
                         << endl;
                    throw "NoN3ToN1Exception";
                }

                double ncx = c6_vec.x(), ncy = c6_vec.y(), ncz = c6_vec.z(), nnx = n1_n3_vec.x(), nny = n1_n3_vec.y(), nnz = n1_n3_vec.z();

                auto z = (ncz*pow(nnx,2) + ncz*pow(nny,2) - ncx*nnx*nnz - ncy*nny*nnz)/
                         (ncx*nnx + ncy*nny);

                z_trans += {0, 0, z / range};

//                auto theta = -acos(-((ncx*nnx + ncy*nny)/
//                                  sqrt(pow(ncx,2)*pow(nnx,2) + pow(ncy,2)*pow(nnx,2) +
//                                       pow(ncx,2)*pow(nny,2) + pow(ncy,2)*pow(nny,2))));
//                double theta = twist;
//                cout << "calculated theta: " << fmod(theta * RAD_TO_DEG, 360) << ", twist: " << fmod(twist * (n + 1) * RAD_TO_DEG, 360) << endl;
            }

            // TODO add these as part of the conditional entrance into if block to avoid unnecessary calculations
            s_trans_prev = hp.getStepTranslationVec(1);
            s_rot_prev = hp.getStepRotationOBMatrix(1);
            g_trans_prev = g_trans;
            g_rot_prev = g_rot;

            //        double twist = 36.5 * DEG_TO_RAD;
            //        double twist = 33.2 * DEG_TO_RAD;
            //        double twist = hp.getTwist() * DEG_TO_RAD;
            //        unsigned n = chain_length_ + 1;
            //        auto z_rot = matrix3x3{{cos(twist * n), -sin(twist * n), 0},{sin(twist * n), cos(twist * n), 0}, {0, 0, 1}};

            // TODO make this more efficient
            // Reverse the tuples back, probably not necessary
            std::reverse(c_n1_or_n3_tuple_.begin(), c_n1_or_n3_tuple_.end());

        }

        // TODO the twist and translation cannot be performed simultaneously, need to do rotation before translation
//        twist = 33.2 * DEG_TO_RAD;
        z_rot = matrix3x3{{cos(twist), -sin(twist), 0},{sin(twist), cos(twist), 0}, {0, 0, 1}};
        for (unsigned i = 0; i < l; ++i) {
            vector3 v3;
            v3.Set(xyz + 3 * i + shift_index);
            v3 += z_trans;
//            v3 *= z_rot;
            v3.Get(xyz + 3 * i + shift_index);
        }
//        for (unsigned i = 0; i < l; ++i) {
//            vector3 v3;
//            v3.Set(xyz + 3 * i + shift_index);
//            v3 *= z_rot;
//            v3.Get(xyz + 3 * i + shift_index);
//        }

    }
}

double Chain::calcRotation(double init_theta, double tol, unsigned max_iterations, vector3 c6_vec, vector3 n1_vec, vector3 n3_vec) {
    double theta = init_theta - 0.1 * DEG_TO_RAD;
    double best_theta = theta;
    double best_angle = 10000;
    auto n1_vec_original = n1_vec;

//    cout << "C6: " << c6_vec << endl;
//    cout << "N1: " << n1_vec_original << endl;
//    cout << "N3: " << n3_vec << endl;

    for (unsigned i = 0; i < max_iterations; ++i) {

        theta = theta + 0.1 * DEG_TO_RAD;
        matrix3x3 m_twist = matrix3x3{{cos(theta), -sin(theta), 0},{sin(theta), cos(theta), 0}, {0, 0, 1}};
        n1_vec = m_twist*n1_vec_original;
        double c6x = c6_vec.x(), c6y = c6_vec.y(), c6z = c6_vec.z(), n1x = n1_vec.x(), n1y = n1_vec.y(), n1z = n1_vec.z(),
                n3x = n3_vec.x(), n3y = n3_vec.y(), n3z = n3_vec.z();


        double angle = RAD_TO_DEG*acos(((n1x - n3x)*(-c6x + n3x) + (n1y - n3y)*(-c6y + n3y))/
                              (sqrt(pow(-c6x + n3x,2) + pow(-c6y + n3y,2))*
                               sqrt(pow(n1x - n3x,2) + pow(n1y - n3y,2))));

//        cout << "theta, angle: " << theta*RAD_TO_DEG << ", " << angle << endl;
//        angle -= 180;
        if (abs(angle) < best_angle) {
            best_theta = theta;
            best_angle = abs(angle);
//            cout << "theta, angle: " << best_theta*RAD_TO_DEG << ", " << best_angle << endl;
        }

//        cout << "Current best angle: " << RAD_TO_DEG*best_angle << endl;

//        cout << "n1_vec: " << n1_vec << endl;
//        cout << "n3_vec: " << n3_vec << endl;
//        cout << "c6_vec: " << c6_vec << endl;

//        double f = acos(((n1z - n3z)*(-c6z + n3z) +
//                           (-c6y + n3y)*(-n3y + n1y*cos(theta) + n1x*sin(theta)) +
//                           (c6x - n3x)*(n3x - n1x*cos(theta) + n1y*sin(theta)))/
//                          (sqrt(pow(c6x - n3x,2) + pow(c6y - n3y,2) + pow(c6z - n3z,2))*
//                           sqrt(pow(n1z - n3z,2) + pow(-n3y + n1y*cos(theta) + n1x*sin(theta),2) +
//                                pow(n3x - n1x*cos(theta) + n1y*sin(theta),2))));
//
//        if (f < tol)
//            return theta;
//
//        double df = -(2*(c6z*n1z - c6x*n3x + pow(n3x,2) - c6y*n3y + pow(n3y,2) - c6z*n3z - n1z*n3z +
//                         pow(n3z,2) + (c6x*n1x + c6y*n1y - n1x*n3x - n1y*n3y)*cos(theta) +
//                         (c6y*n1x - c6x*n1y + n1y*n3x - n1x*n3y)*sin(theta))*
//                      ((n1y*n3x - n1x*n3y)*cos(theta) + (n1x*n3x + n1y*n3y)*sin(theta)) +
//                      2*((c6x - n3x)*(n1y*cos(theta) + n1x*sin(theta)) +
//                         (-c6y + n3y)*(n1x*cos(theta) - n1y*sin(theta)))*
//                      (pow(n1z - n3z,2) + pow(-n3y + n1y*cos(theta) + n1x*sin(theta),2) +
//                       pow(n3x - n1x*cos(theta) + n1y*sin(theta),2)))/
//                    (2.*sqrt(pow(c6x - n3x,2) + pow(c6y - n3y,2) + pow(c6z - n3z,2))*
//                     pow(pow(n1z - n3z,2) + pow(-n3y + n1y*cos(theta) + n1x*sin(theta),2) +
//                           pow(n3x - n1x*cos(theta) + n1y*sin(theta),2),1.5)*
//                     sqrt(1 - pow((n1z - n3z)*(-c6z + n3z) +
//                                    (-c6y + n3y)*(-n3y + n1y*cos(theta) + n1x*sin(theta)) +
//                                    (c6x - n3x)*(n3x - n1x*cos(theta) + n1y*sin(theta)),2)/
//                              ((pow(c6x - n3x,2) + pow(c6y - n3y,2) + pow(c6z - n3z,2))*
//                               (pow(n1z - n3z,2) + pow(-n3y + n1y*cos(theta) + n1x*sin(theta),2) +
//                                pow(n3x - n1x*cos(theta) + n1y*sin(theta),2)))));
//
//        theta = theta - f / df;
    }
//    cout << "Best theta, best angle: " << best_theta*RAD_TO_DEG << ", " << best_angle << endl;
    return best_theta;
//    cerr << "Failed to converge" << endl;
//    return 0;
}
