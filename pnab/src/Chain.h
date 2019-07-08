//
// Created by jbarnett8 on 2/10/18.
//

#ifndef PNAB_CHAIN_H
#define PNAB_CHAIN_H

#include <openbabel/forcefield.h>
#include "Containers.h"

#define KJ_TO_KCAL 0.239006

class Chain {

public:
    Chain() {};
    Chain(PNAB::Bases bases, PNAB::Backbone backbone, std::vector<std::string> strand,
          std::string ff_type, std::array<unsigned, 2> &range, bool double_stranded = true, bool hexad = false);
    ~Chain() {
        for (auto v : base_coords_vec_)
            delete v;
    }

    PNAB::ConformerData generateConformerData(double *xyz, PNAB::HelicalParameters &hp);

    OpenBabel::OBMol getChain() {
        if (double_stranded_ || hexad_)
            return combined_chain_;
        else
            return chain_;
    }

private:
    OpenBabel::OBMol chain_, c_chain_, d_chain_, e_chain_, f_chain_, g_chain_, combined_chain_;
    std::vector<unsigned> new_bond_ids_, c_new_bond_ids_, d_new_bond_ids_, e_new_bond_ids_, f_new_bond_ids_, g_new_bond_ids_;
    std::vector<unsigned> deleted_atoms_ids_, c_deleted_atoms_id_, d_deleted_atoms_id_, e_deleted_atoms_id_, f_deleted_atoms_id_, g_deleted_atoms_id_;
    std::vector<unsigned> num_bu_A_mol_atoms_, num_bu_B_mol_atoms_, num_bu_C_mol_atoms_, num_bu_D_mol_atoms, num_bu_E_mol_atoms, num_bu_F_mol_atoms;
    std::vector<std::size_t> base_connect_index, c_base_connect_index, d_base_connect_index, e_base_connect_index, f_base_connect_index, g_base_connect_index;
    std::vector<double*> base_coords_vec_, c_base_coords_vec_, d_base_coords_vec_, e_base_coords_vec_, f_base_coords_vec_, g_base_coords_vec_;
    unsigned chain_length_;
    bool isKCAL_, double_stranded_, hexad_;
    OpenBabel::OBForceField *pFF_;
    std::array<unsigned, 2> monomer_bb_index_range_;
    std::vector<unsigned> bb_start_index_, c_bb_start_index_, d_bb_start_index_, e_bb_start_index_, f_bb_start_index_, g_bb_start_index_;
    std::string ff_type_;
    PNAB::Backbone backbone_;
    std::vector<std::tuple<std::string,unsigned,unsigned>> n1_or_n3_tuple_, c_n1_or_n3_tuple_;
    OpenBabel::vector3 z_trans;

    OpenBabel::OBFFConstraints constraintsTot_, constraintsTor_, constraintsAng_, constraintsBond_;

    void fillConformerEnergyData(double *xyz, PNAB::ConformerData &conf_data);
    void setupChain(std::vector<PNAB::Base> &strand, OpenBabel::OBMol &chain, std::vector<unsigned> &new_bond_ids,
                        std::vector<unsigned> &deleted_atoms_ids, std::vector<unsigned> &num_base_unit_atoms,
                        std::vector<unsigned> &bb_start_index, std::vector<double *> &base_coords_vec,
                        unsigned chain_index, bool double_stranded, bool hexad);
    void setupFFConstraints(OpenBabel::OBMol &chain, std::vector<unsigned> &new_bond_ids, unsigned offset = 0);
    void setCoordsForChain(double *xyz, double *conf, PNAB::HelicalParameters &hp, std::vector<unsigned> &num_bu_atoms,
                           std::vector<unsigned> &bb_start_index, std::vector<double *> &base_coords_vec,
                           std::vector<unsigned> &deleted_atoms_ids, bool is_second_strand, bool is_hexad, unsigned chain_index);

};


#endif //PNAB_CHAIN_H
