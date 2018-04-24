//
// Created by jbarnett8 on 2/26/18.
//

#ifndef PNAB_MONTECARLOROTORSEARCH_H
#define PNAB_MONTECARLOROTORSEARCH_H

#include <openbabel/obconversion.h>
#include <openbabel/forcefield.h>
#include "../FileParser.h"
#include "../UnitChain.h"
#include "../Chain.h"

class MonteCarloRotorSearch {

public:
    MonteCarloRotorSearch(PNAB::RuntimeParameters &runtime_params, PNAB::Base base_a, PNAB::Base base_b,
                          PNAB::Backbone backbone, PNAB::HelicalParameters &helical_params);

    bool run();
    OpenBabel::OBMol getMolecule() {
        return mol;
    }
    std::vector<double*> getCoordinates() {
        return coord_vec;
    }

private:
    PNAB::RuntimeParameters runtime_params_;
    PNAB::Base base_a_, base_b_;
    PNAB::Backbone backbone_;
    PNAB::HelicalParameters helical_params_;
    std::mt19937_64 rng;
    std::vector<PNAB::ConformerData> conf_data_vec;
    std::vector<double *> coord_vec;
    OpenBabel::matrix3x3 step_rot;
    std::array<double, 9> glbl_rot;
    OpenBabel::vector3 step_translate, glbl_translate;
    OpenBabel::OBMol mol;
    OpenBabel::OBConversion conv;
    OpenBabel::OBMol test_chain;
    unsigned monomer_num_coords;

    double measureDistance(double *coords, unsigned head, unsigned tail);
    bool isPassingEFilter(const PNAB::ConformerData &conf_data);
    void print(PNAB::ConformerData conf_data);
    double calcRMSD(double *ref, double *conf, unsigned long size) {
        double rmsd = 0;
        for (int i = 0; i < size; i++) {
            rmsd += pow(*(ref + i) - *(conf + i),2);
        }
        return sqrt(rmsd/size);
    }
};

#endif //PNAB_MONTECARLOROTORSEARCH_H
