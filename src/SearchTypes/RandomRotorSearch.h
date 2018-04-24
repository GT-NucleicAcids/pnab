//
// Created by jbarnett8 on 9/17/17.
//

#ifndef PNAB_RANDOMROTORSEARCH_H
#define PNAB_RANDOMROTORSEARCH_H

#include "../FileParser.h"
#include "../UnitChain.h"

#include <random>
#include <functional>

/**
 * \brief A small class to hold a vector of functions. Idea shamelessly taken from
 * https://stackoverflow.com/questions/15128444/c-calling-a-function-from-a-vector-of-function-pointers-inside-a-class-where-t
 */
class VectorOfFunctions {
public:
    template<typename Function>
    void push_back(Function && f) {
        functions.push_back(std::forward<Function>(f));
    }

    void invoke_all() {
        for(auto && f : functions)
            f();
    }


private:
    std::vector<std::function<void()>> functions;
};

class RandomRotorSearch {
public:
    RandomRotorSearch(PNAB::RuntimeParameters &runtime_params, PNAB::UnitChain &uc);
    OpenBabel::OBMol search();
    void speak() {
        std::cout << "Random Rotor Search Speaks!" << std::endl;
    }

    bool passesBackboneLinkerDistanceFilter(double *coord);

    void updateGeometry();

private:
    PNAB::RuntimeParameters runtime_params;
    PNAB::UnitChain uc;
    OpenBabel::OBMol chain;
    std::vector<std::array<std::size_t, 2>> backbone_interlinks;
    double max_backbone_interlink_distance;
    std::array< double, 3 > rotation, translation;
    std::mt19937_64 rng;
    std::vector<std::uniform_real_distribution<double>> dist_vec;
    //std::vector<PNAB::RuntimeParameters::RangeField> range_vec;
    OpenBabel::matrix3x3 uc_twist;
    OpenBabel::vector3 uc_rise;
    std::vector<double*> coords_vec;

    std::vector<std::function<void()>> range_params_functions;
    VectorOfFunctions params_functions;


    enum RangeTypes {
        X_TRANSLATE = 0,
        Y_TRANSLATE,
        Z_TRANSLATE,
        X_ROTATION,
        Y_ROTATION,
        Z_ROTATION
    };


    //void addRange(PNAB::RuntimeParameters::RangeField &field, double &val);
};

#endif //PNAB_RANDOMROTORSEARCH_H
