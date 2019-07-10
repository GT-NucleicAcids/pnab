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
    Chain(PNAB::Bases bases, const PNAB::Backbone &backbone, std::vector<std::string> strand,
          std::string ff_type, std::array<unsigned, 2> &range, bool double_stranded = true, bool hexad = false);
    ~Chain() {
        for (auto v : v_base_coords_vec_) {
            for (auto i : v) 
                delete i;
        }
    }

    PNAB::ConformerData generateConformerData(double *xyz, PNAB::HelicalParameters &hp);

    OpenBabel::OBMol getChain() {
        if (double_stranded_ || hexad_)
            return combined_chain_;
        else
            return v_chain_[0];
    }

private:
    std::vector<OpenBabel::OBMol> v_chain_ = std::vector<OpenBabel::OBMol>(6); //declare 6 chains in case hexad is requested
    OpenBabel::OBMol combined_chain_;
    std::vector<std::vector<unsigned>> v_new_bond_ids_ = std::vector<std::vector<unsigned>>(6);
    std::vector<std::vector<unsigned>> v_deleted_atoms_ids_ = std::vector<std::vector<unsigned>>(6);
    std::vector<std::vector<unsigned>> v_num_bu_A_mol_atoms_ = std::vector<std::vector<unsigned>>(6);
    std::vector<std::vector<std::size_t>> v_base_connect_index = std::vector<std::vector<size_t>>(6);
    std::vector<std::vector<double*>> v_base_coords_vec_ = std::vector<std::vector<double*>>(6);
    unsigned chain_length_, n_chains_;
    bool isKCAL_, double_stranded_, hexad_;
    OpenBabel::OBForceField *pFF_;
    std::array<unsigned, 2> monomer_bb_index_range_;
    std::vector<std::vector<unsigned>> v_bb_start_index_ = std::vector<std::vector<unsigned>>(6);
    std::string ff_type_;
    PNAB::Backbone backbone_;

    OpenBabel::OBFFConstraints constraintsTot_, constraintsTor_, constraintsAng_, constraintsBond_;

    void fillConformerEnergyData(double *xyz, PNAB::ConformerData &conf_data);
    void setupChain(std::vector<PNAB::Base> &strand, OpenBabel::OBMol &chain, std::vector<unsigned> &new_bond_ids,
                        std::vector<unsigned> &deleted_atoms_ids, std::vector<unsigned> &num_base_unit_atoms,
                        std::vector<unsigned> &bb_start_index, std::vector<double *> &base_coords_vec,
                        unsigned chain_index);
    void setupFFConstraints(OpenBabel::OBMol &chain, std::vector<unsigned> &new_bond_ids, unsigned offset = 0);
    void setCoordsForChain(double *xyz, double *conf, PNAB::HelicalParameters &hp, std::vector<unsigned> &num_bu_atoms,
                           std::vector<unsigned> &bb_start_index, std::vector<double *> &base_coords_vec,
                           std::vector<unsigned> &deleted_atoms_ids, unsigned chain_index);

};


#endif //PNAB_CHAIN_H
