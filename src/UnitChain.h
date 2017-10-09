//
// Created by jbarnett8 on 9/4/17.
//

#ifndef PNAB_UNITCHAIN_H
#define PNAB_UNITCHAIN_H

#include "Containers.h"

namespace PNAB {

/**
 * The single most basic unit used for testing conformations. Scenarios where periodic boundary conditions
 * are satisfied can have as little as a single \code{Base} in the chain while sets of bases that do not
 * satisfy this condition may require more. PNAB always assumes periodic boundary conditions between
 * \code{UnitChain}s.
 *
 * A unit chain is constructed from a list of strings containing codes matching to provided bases. The
 * length of the \code{chain_string} vector is the symmetry of the \code{UnitChain} such that the first
 * element of each string makes the first "rung" in the chain.
 */
    class UnitChain {

    private:
        std::vector< std::vector< PNAB::BaseUnit > > base_units;
        OpenBabel::OBMol backbone;

    public:
        UnitChain(std::vector<std::string> chain_string, Bases bases, PNAB::Backbone backbone);
        PNAB::BaseUnit getUnit(std::size_t strand, std::size_t index) {
            return base_units[strand][index];
        }
    };

}

#endif //PNAB_UNITCHAIN_H
