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
    Chain(PNAB::Bases bases, PNAB::Backbone backbone, std::vector<std::string> strand,
          std::string ff_type, std::array<unsigned, 2> &range, bool double_stranded = true);
    ~Chain() {
        delete constraintsAng_, constraintsBond_, constraintsTor_, constraintsTot_;
        for (auto v : base_coords_vec_)
            delete v;
    }

    PNAB::ConformerData generateConformerData(double *xyz, PNAB::HelicalParameters &hp);

    OpenBabel::OBMol getChain() {
        return chain_;
    }

private:
    OpenBabel::OBMol chain_;
    std::vector<unsigned> newBondIDs_;
    std::vector<unsigned> deletedAtomsId_;
    std::vector<unsigned> num_bu_A_mol_atoms_;
    std::vector<double*> base_coords_vec_;
    unsigned chain_length_;
    bool isKCAL_, double_stranded_;
    OpenBabel::OBForceField *pFF_;
    std::array<unsigned, 2> monomer_bb_index_range_;
    std::vector<unsigned> bb_start_index_;
    std::string ff_type_;

    OpenBabel::OBFFConstraints *constraintsTot_, *constraintsTor_, *constraintsAng_, *constraintsBond_;

    void fillConformerEnergyData(double *xyz, PNAB::ConformerData &conf_data);

};


#endif //PNAB_CHAIN_H
