//
// Created by jbarnett8 on 2/26/18.
//

#ifndef PNAB_MONTECARLOROTORSEARCH_H
#define PNAB_MONTECARLOROTORSEARCH_H

#include <openbabel/obconversion.h>
#include <openbabel/forcefield.h>
#include "../FileParser.h"
#include "../Chain.h"
#include <openbabel/rotor.h>

class MonteCarloRotorSearch {

public:
    /**
     * \brief High level rotor search function used to find ideal conformation of arbitrary backbone and helical
     * parameter combinations
     * @param runtime_params Series of parameters that controls how the search algorithm runs
     * @param backbone The molecule over which the algorithm searches
     * @param helical_params The geometric parameters constraining possible conformations of the backbone
     * @param bases List of possible bases that serves of as the library used by 'strand'
     * @param strand List of base names that gives the identity of the test chain
     */
    MonteCarloRotorSearch(PNAB::RuntimeParameters &runtime_params, PNAB::Backbone backbone,
                          PNAB::HelicalParameters &helical_params, PNAB::Bases bases, std::vector<std::string> strand);

    bool run();

private:
    PNAB::RuntimeParameters runtime_params_;
    PNAB::Base base_a_;
    std::array<unsigned, 2> backbone_range_;
    PNAB::Backbone backbone_;
    PNAB::HelicalParameters helical_params_;
    PNAB::Bases bases_;
    std::vector<std::string> strand_;
    std::mt19937_64 rng_;
    std::vector<PNAB::ConformerData> conf_data_vec_;
    OpenBabel::matrix3x3 step_rot_;
    std::array<double, 9> glbl_rot_;
    OpenBabel::vector3 step_translate_, glbl_translate_;
    OpenBabel::OBConversion conv_;
    OpenBabel::OBMol test_chain_;
    unsigned monomer_num_coords_;

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
