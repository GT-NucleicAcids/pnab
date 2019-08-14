//
// Created by jbarnett8 on 2/26/18.
//

#ifndef PNAB_CONFORMATIONSEARCH_H
#define PNAB_CONFORMATIONSEARCH_H

#include <openbabel/obconversion.h>
#include <openbabel/forcefield.h>
#include "../Chain.h"
#include <openbabel/rotor.h>

class ConformationSearch {

public:
    /**
     * \brief High level rotor search function used to find ideal conformation of arbitrary backbone and helical
     * parameter combinations
     * @param runtime_params Series of parameters that controls how the search algorithm runs
     * @param backbone The molecule over which the algorithm searches
     * @param helical_params The geometric parameters constraining possible conformations of the backbone
     * @param bases List of possible bases that serves of as the library used by 'strand'
     * @param strand List of base names that gives the identity of the test chain
     * @param prefix to conformer index for multiprocess computations
     */
    ConformationSearch(PNAB::RuntimeParameters &runtime_params, PNAB::Backbone backbone,
                       PNAB::HelicalParameters &helical_params, PNAB::Bases bases,
                       std::string prefix);

    std::string run();
    std::string RandomSearch(bool weighted);
    std::string MonteCarloSearch(bool weighted);
    std::string SystematicSearch();
    std::vector <std::piecewise_linear_distribution<double>> weighted_distributions(OpenBabel::OBRotorList &rl, OpenBabel::OBMol &bu_a_mol);

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
    OpenBabel::matrix3x3 step_rot_, glbl_rot_;
    OpenBabel::vector3 step_translate_, glbl_translate_;
    OpenBabel::OBConversion conv_;
    OpenBabel::OBMol test_chain_;
    unsigned monomer_num_coords_;
    bool is_double_stranded_, is_hexad_, is_parallel_;
    std::string ff_type_;
    std::string prefix_;

    double measureDistance(double *coords, unsigned head, unsigned tail);
    std::string print(PNAB::ConformerData conf_data);
    double calcRMSD(double *ref, double *conf, unsigned long size) {
        double rmsd = 0;
        for (int i = 0; i < size; i++) {
            rmsd += pow(*(ref + i) - *(conf + i),2);
        }
        return sqrt(rmsd/size);
    }
    void printProgress(std::size_t search_index, std::size_t search_size);
};

#endif //PNAB_MONTECARLOROTORSEARCH_H
