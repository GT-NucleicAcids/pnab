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
          std::string ff_type, std::array<unsigned, 2> &range, bool double_stranded = true);
    ~Chain() {
        for (auto v : base_coords_vec_)
            delete v;
    }

    PNAB::ConformerData generateConformerData(double *xyz, PNAB::HelicalParameters &hp);

    OpenBabel::OBMol getChain() {
        if (double_stranded_)
            return combined_chain_;
        else
            return chain_;
    }

private:
    OpenBabel::OBMol chain_, c_chain_, combined_chain_;
    std::vector<unsigned> new_bond_ids_, c_new_bond_ids_;
    std::vector<unsigned> deleted_atoms_ids_, c_deleted_atoms_id_;
    std::vector<unsigned> num_bu_A_mol_atoms_, num_bu_B_mol_atoms_;
    std::vector<std::size_t> base_connect_index, c_base_connect_index;
    std::vector<double*> base_coords_vec_, c_base_coords_vec_;
    unsigned chain_length_;
    bool isKCAL_, double_stranded_;
    OpenBabel::OBForceField *pFF_;
    std::array<unsigned, 2> monomer_bb_index_range_;
    std::vector<unsigned> bb_start_index_, c_bb_start_index_;
    std::string ff_type_;
    PNAB::Backbone backbone_;
    std::vector<std::tuple<std::string,unsigned,unsigned>> n1_or_n3_tuple_, c_n1_or_n3_tuple_;
    OpenBabel::vector3 z_trans;

    OpenBabel::OBFFConstraints constraintsTot_, constraintsTor_, constraintsAng_, constraintsBond_;

    void fillConformerEnergyData(double *xyz, PNAB::ConformerData &conf_data);
    void setupChain(std::vector<PNAB::Base> &strand, OpenBabel::OBMol &chain, std::vector<unsigned> &new_bond_ids,
                        std::vector<unsigned> &deleted_atoms_ids, std::vector<unsigned> &num_base_unit_atoms,
                        std::vector<unsigned> &bb_start_index, std::vector<double *> &base_coords_vec,
                        std::vector<std::tuple<std::string, unsigned, unsigned>> &h1_or_h3_tuple,
                        unsigned chain_index, bool double_stranded);
    void setupFFConstraints(OpenBabel::OBMol &chain, std::vector<unsigned> &new_bond_ids, unsigned offset = 0);
    void setCoordsForChain(double *xyz, double *conf, PNAB::HelicalParameters &hp, std::vector<unsigned> &num_bu_atoms,
                           std::vector<unsigned> &bb_start_index, std::vector<double *> &base_coords_vec,
                           std::vector<unsigned> &deleted_atoms_ids, bool is_second_strand);

    double calcRotation(double init_z, double tol, unsigned max_iterations,OpenBabel::vector3 c6_vec,
                        OpenBabel::vector3 n1_vec, OpenBabel::vector3 n3_vec);

};


#endif //PNAB_CHAIN_H
