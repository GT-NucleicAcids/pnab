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
    Chain(PNAB::BaseUnit bu, unsigned chain_length, std::string ff_type);
    ~Chain() {
        delete constraintsAng, constraintsBond, constraintsTor, constraintsTot;
    }

    PNAB::ConformerData generateConformerData(double *xyz, PNAB::HelicalParameters hp);

    OpenBabel::OBMol getChain() {
        return chain;
    }

private:
    OpenBabel::OBMol monomer;
    OpenBabel::OBMol chain;
    std::vector<unsigned> newBondIDs;
    std::vector<unsigned> deletedAtomsId;
    unsigned chain_length_;
    bool isKCAL;
    OpenBabel::OBForceField *pFF;

    OpenBabel::OBFFConstraints *constraintsTot, *constraintsTor, *constraintsAng, *constraintsBond;

    void fillConformerEnergyData(double *xyz, PNAB::ConformerData &conf_data);

};


#endif //PNAB_CHAIN_H
