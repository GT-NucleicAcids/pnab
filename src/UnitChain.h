//
// Created by jbarnett8 on 9/4/17.
//

#ifndef PNAB_UNITCHAIN_H
#define PNAB_UNITCHAIN_H

#include "Containers.h"
#include <openbabel/rotor.h>

namespace PNAB {

/**
 * \brief The single most basic unit used for testing conformations.
 *
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
        /**
         * \brief Adds the base and backbone monomer to the longer chain defining the UnitChain
         * @param bu The BaseUnit to add to the chain member
         */
        void addToChain(BaseUnit &bu);

        std::vector< std::vector< PNAB::BaseUnit > > base_units;            //!< \brief 2D vector of BaseUnit values organized as [strand][index along strand]
        std::vector< std::array<std::size_t, 2 > > base_indices,            //!< \brief Indices indicating ranges that are bases in the form {[start, stop], [start, stop], ...}
                                                   backbone_indices,        //!< \brief Indices indicating ranges that are backbones in the form {[start, stop], [start, stop], ...}
                                                   base_unit_indices,       //!< \brief Indices indicating start and stop of specific BaseUnit element in form {[start, stop], [start, stop], ...}
                                                   backbone_linker_indices, //!< \brief Indices indicating the linkers between successive backbones
                                                   base_bond_indices;       //!< \brief Indices for start and stop of bond indices in each base
        OpenBabel::OBMol backbone, chain;                                   //!< \brief Raw construction of chain from BaseUnit elements in base_units
        std::vector<double> original_coords;                                //!< \brief Used to store original coordinates prior to geometric modifications so we can reset at any point

    public:

        UnitChain() {}

        /**
         * Generates a member of the UnitChain class
         * @param chain_string A vector of strings that contain the three letter code matching in the list of Base members, bases
         * @param bases A vector of available Base members
         * @param backbone The backbone used to attached to the bases
         */
        UnitChain(std::vector<std::string> chain_string, Bases bases, PNAB::Backbone backbone);
        PNAB::BaseUnit getUnit(std::size_t strand, std::size_t index) {
            return base_units[strand][index];
        }

        const std::vector< std::array<std::size_t, 2 > >& getBaseIndices() {
            return base_indices;
        };

        const std::vector< std::array<std::size_t, 2 > >& getBackboneLinkerIndices() {
            return backbone_linker_indices;
        };

        const std::vector< std::array<std::size_t, 2 > >& getBaseBondIndices() {
            return base_bond_indices;
        };

        double* getCoordinates() {
            return chain.GetCoordinates();
        }

        /**
         * \brief Update coordinates of raw chain based off of geometric parameters provided
         * @param translation How BaseUnit elements are to be translated in the form {x_displacement, y_displacement, rise} in Angstroms
         * @param rotation How BaseUnit elements are to be rotated in the form {inclination, tip, twist} in degrees
         */
        void updateHelicalParameters(std::array< double, 3 > &translation, std::array< double, 3 > &rotation);

        /**
         * Gives a copy of the current raw chain OBMol object
         * @return The raw chain
         */
        OpenBabel::OBMol getChainMol() {
            return chain;
        }
    };

}

#endif //PNAB_UNITCHAIN_H
