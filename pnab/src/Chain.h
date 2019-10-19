/**@file
 * @brief A file for declaring a class for building and evaluating nucleic acid strands
 */

#ifndef PNAB_CHAIN_H
#define PNAB_CHAIN_H

#include <algorithm>
#include <openbabel/forcefield.h>
#include <openbabel/rotor.h>
#include <openbabel/generic.h>
#include <openbabel/obiter.h>
#include "Containers.h"
#define KJ_TO_KCAL 0.239006

namespace PNAB {

    /**
    * @brief A class for building nucleic acid strands and evaluating their energies
    *
    * The class creates the strands by connecting the nucleotides created in the BaseUnit
    * class. It forms bonds between the nucleotides for the given sequence of the nucleobases,
    * and creates duplex or hexad systems if requested. The class also contains functions
    * for computing the energy terms for the system and checking whether the candidates are
    * accepted. Several functions in the class identifies the proper energy terms to be computed
    * for the system (e.g. bond and angle energies).
    *
    * @sa ConformationSearch
    * @sa HelicalParameters
    * @sa RuntimeParameters
    * @sa BaseUnit
    * @sa Bases
    */
    class Chain {

    public:
        /**
        * @brief Constructor for the chain class
        *
        * The constructor performs several tasks. It creates an OpenBabel::ForceField object and check
        * that it works as expected. It identifies whether we are working with the canonical nucleobases
        * or with the hexad. For double-stranded and hexad systems, it checks whether the user provided
        * correct comlimentary bases. After all these checks, it proceeds to setup the requested molecules.
        * The constructor does not put correct coordinates for the molecules, but creates the topology
        * of each strand in the system. The constructor also identifies the atoms that should be ignored when
        * computing the separate energies (e.g. bond energy) for the backbone candidates.
        *
        * @param bases A vector of the all the defined bases,
        * @param backbone The backbone
        * @param strand A vector of the names of all bases in the strand, RuntimeParameters::strand
        * @param ff_type The force field type, RuntimeParameters::ff_type
        * @param range Backbone index range for the first nucleotide
        * @param hexad Defines whether the 60 degrees rotation for hexads is performed, RuntimeParameters::is_hexad
        * @param build_strand Defines whether to build a given strand, RuntimeParameters::build_strand
        * @param strand_orientation The orientation of each strand in the hexad, RuntimeParameters::strand_orientation
        *
        * @sa setupChain
        * @sa setupFFConstraints
        */
        Chain(PNAB::Bases bases, const PNAB::Backbone &backbone, std::vector<std::string> strand,
              std::string ff_type, std::array<unsigned, 2> &range, bool hexad,
              std::vector<bool> build_strand = {true, false, false, false, false, false},
              std::vector<bool> strand_orientation = {true, true, true, true, true, true});

        /**
        * @brief Destructor for the chain class
        */
        ~Chain() {
            for (auto v : v_base_coords_vec_) {
                for (auto i : v)
                    delete i;
            }
        }

        /**
        * @brief Generate structure and energy data for nucleic acid conformers
        *
        * After the conformation search finds a candidate that satisfies the distance threshold,
        * this function is called with the coordinates of the backbone candidate to build the full system
        * and computes its energy properties.
        *
        * This function calls Chain::setCoordsForChain and Chain::fillConformerEnergyData to setup the coordinates and
        * compute the energies, respectively.
        *
        * @param xyz The coordinates of the backbone candidate that satisfied the distance threshold
        * @param hp An instance of HelicalParameters that contains the helical structure data and functions
        * @param energy_filter The energy thresholds for the system, RuntimeParameters::energy_filter
        *
        * @returns The conformer energy data, ConformerData, containing whether the candidate is accepted
        *
        * @sa setCoordsForChain
        * @sa fillConformerEnergyData
        * @sa ConformationSearch::RandomSearch
        * @sa ConformationSearch::MonteCarloSearch
        * @sa ConformationSearch::SystematicSearch
        * @sa ConformationSearch::GeneticAlgorithmSearch
        */
        PNAB::ConformerData generateConformerData(double *xyz, PNAB::HelicalParameters &hp, std::vector<double> energy_filter);

        /**
         * @brief Returns the openbabel molecule containing the structure of the system
         *
         * @returns The openbabel molecule containing the structure of the system
         */
        OpenBabel::OBMol getChain() {
            return combined_chain_;
        }

    private:
        std::vector<OpenBabel::OBMol> v_chain_ = std::vector<OpenBabel::OBMol>(6); //!< @brief A vector of OpenBabel::OBMol containing the molecules for each strand in the system
        OpenBabel::OBMol combined_chain_; //!< @brief An OpenBabel::OBMol molecule containing the structure of the whole system
        std::vector<std::vector<unsigned>> v_new_bond_ids_ = std::vector<std::vector<unsigned>>(6); /*!< @brief a vector containing a vector of the IDs of the atoms forming
                                                                                                     *    new bonds between the nucleotides in each strand
                                                                                                     */
        std::vector<std::vector<unsigned>> v_deleted_atoms_ids_ = std::vector<std::vector<unsigned>>(6); /*!< @brief A vector containing a vector of the IDs of the atoms
                                                                                                          *   deleted in each strand because of the formation of new bonds
                                                                                                          */
        std::vector<std::vector<unsigned>> v_num_bu_A_mol_atoms_ = std::vector<std::vector<unsigned>>(6); /*!< @brief A vector containing a vector of the number of atoms
                                                                                                           *   in each BaseUnit for each strand
                                                                                                           */
        std::vector<std::vector<std::vector<unsigned>>> v_fixed_bonds = std::vector<std::vector<std::vector<unsigned>>>(6); /*!< @brief A vector containing the indices of
                                                                                                                             *   fixed rotatable bonds for each strad
                                                                                                                             */
        std::vector<std::vector<double*>> v_base_coords_vec_ = std::vector<std::vector<double*>>(6); /*!< @brief A vector containing a vector the coordinates of
                                                                                                      *   each nucleotide in each strand
                                                                                                      */
        unsigned chain_length_, //!< @brief The number of nucleotides in the strand
                 n_chains_; //!< @brief The number of strands in the system
        bool isKCAL_, //!< @brief Whether the energy computed by openbabel is in kcal/mol
             hexad_; //!< @brief Whether we are building a hexad, RuntimeParameters::is_hexad
        std::vector<bool> strand_orientation_; //!< @brief A vector containing the orientation of each strand in the hexad, RuntimeParameters::strand_orientation
        OpenBabel::OBForceField *pFF_; //!< @brief The openbabel force field. Used to compute the energy of the system
        std::array<unsigned, 2> monomer_bb_index_range_; //!< @brief Backbone index range for the first nucleotide
        std::vector<std::vector<unsigned>> v_bb_start_index_ = std::vector<std::vector<unsigned>>(6); /*!< @brief A vector containing a vector of the starting
                                                                                                       *   indices of the backbone atoms in the BaseUnit for each strand
                                                                                                       */
        std::string ff_type_; //!< @brief The force field type (e.g. "GAFF"), RuntimeParameters::ff_type
        std::vector<std::vector<unsigned int>> all_angles_, //!< @brief A vector of the vector of atom indices forming all the angles for which we need to compute the energy
                                               all_torsions_; //!< @brief A vector of the vector of atom indices forming all the torsions for which we need to compute the energy
        std::vector<bool> is_fixed_bond; //!< @brief A vector containing whether the torsional energy term is for a fixed rotatable bond or not
        std::vector<bool> build_strand_; //!< @brief A vector containing whether a given strand should be built

        OpenBabel::OBFFConstraints constraintsBond_, //!< @brief Setting all atoms not forming the new bond between the first two nucleotides to be ignored during bond energy computation
                                   constraintsAng_, //!< @brief An empty constraint object for angles; Energy groups are used for the angle terms
                                   constraintsTor_, //!< @brief An empty constraint object for torsions; Energy groups are used for the torsion terms
                                   constraintsTot_; //!< @brief An empty constraint object for van der Waals and total energy terms. No ignored atoms in these computations

        /**
        * @brief Computes the energy terms for the candidate system and determines whether it satisfies the energy thresholds
        *
        * This function uses openbabel to compute the energy for the candidate systems. See the description in setupFFConstraints
        * for the bond, angle, and torsional energies terms. The van der Waals and total energy terms of the systems are computed
        * without any constraints. The energy terms are computed sequentially as follows: bond, angle, torsion, van der Waals, and
        * total energies. If any one of the energy terms does not satisfy the energy thresholds defined in RuntimeParameters::energy_filter,
        * then the candidate is rejected, and we do not proceed to compute the additional energy terms.
        *
        * @param xyz The coordinates of the whole system
        * @param conf_data An object that will hold the values of the energy terms, the coordinates, and whether it is accepted.
        * @param energy_filter The energy thresholds for the system, RuntimeParameters::energy_filter
        *
        * @sa generateConformerData
        */
        void fillConformerEnergyData(double *xyz, PNAB::ConformerData &conf_data, std::vector<double> energy_filter);

        /**
        * @brief Creates the molecule for each strand in the system
        *
        * This function creates the topology (makes and removes bonds) for each strand in the system.
        * It also determines the indices for all the fixed rotatable bonds in the system.
        *
        * @param strand The list of bases in the strand
        * @param chain An empty openbabel molecule to be populated with the topology of the strand
        * @param new_bond_ids An empty vector to be populated with the indices of the atoms forming bonds between the nucleotides
        * @param deleted_atoms_ids An empty vector to be populated with the indices of the atoms deleted because of the bonding between the nucleotides
        * @param num_base_unit_atoms An empty vector to be populated with the number of atoms in each BaseUnit
        * @param bb_start_index An empty vector to be populated with the starting index of the backbone atoms in the BaseUnit
        * @param base_coords_vec An empty vector to be populated with the coordinates of each one of the nucleotides
        * @param fixed_bonds_vec An empty vector to be populated with the indices of fixed bonds in each one of the nucleotides
        * @param backbone An instance of the Backbone class
        * @param chain_index The index of the chain
        *
        * @sa Chain::Chain
        */
        void setupChain(std::vector<PNAB::Base> &strand, OpenBabel::OBMol &chain, std::vector<unsigned> &new_bond_ids,
                            std::vector<unsigned> &deleted_atoms_ids, std::vector<unsigned> &num_base_unit_atoms,
                            std::vector<unsigned> &bb_start_index, std::vector<double *> &base_coords_vec,
                            std::vector<std::vector<unsigned>> &fixed_bonds_vec, const Backbone &backbone, unsigned chain_index);
        /**
        * @brief Determines the terms that should be ignored during the computation of the bond, angle, and torsional energies
        *
        * For the bond energy, we consider only the one new bond that is formed between the first and second nucleobases in the system.
        * This bond energy is the same for all new bonds, so we do not need to compute it for all of them.
        *
        * For the angle energy, we consider all the new angles that are formed between the first and second nucleobases in the system.
        * This angle energy is the same for all new bonds, so we do not need to compute it for all of them.
        *
        * For the torsional energy, we consider all the rotatable bonds in all the nucleobases. There will be some redundant calculations,
        * as some torsional angles will be the same. This code also identifies whether a rotatable torsional angle is fixed or not.
        * For the identified torsional angles, we do not include any torsional angle that runs between the residues, as these are
        * not rotated in the search procedure.
        *
        * @param chain The openbabel molecule that contain the topology of the strand
        * @param new_bond_ids A vector containing the indices of the atoms forming bonds between the nucleotides
        * @param fixed_bonds_vec A vector containing the indices of fixed bonds in each one of the nucleotides
        * @param offset The number of atoms that have been processed previously. Used to set correct indices for all the strands.
        *
        * @sa Chain::Chain
        * @sa setupChain
        * @sa fillConformerEnergyData
        *
        */
        void setupFFConstraints(OpenBabel::OBMol &chain, std::vector<unsigned> &new_bond_ids, std::vector<std::vector<unsigned>> &fixed_bonds_vec, unsigned offset = 0);

        /**
        * @brief Set the coordinates for each strand in the system
        *
        * This function set the coordinates for each nucleotide in the strand. For DNA- and RNA- like systems, this function can
        * generate anti-parallel strands. For hexad systems, it can generate any combination of the parallel and anti-parallel
        * strands.
        *
        * @param xyz An array that will get the coordinates of the strands
        * @param conf An array containing the coordinates of one backbone determined by the search algorithms
        * @param hp An instance of HelicalParameters which has the helical parameters and functions to generate the coordinates
        * @param num_bu_atoms The number of atoms in each BaseUnit
        * @param bb_start_index A vector containing the starting indices of the backbone atoms in the BaseUnit
        * @param base_coords_vec A vector containing the coordinates of each one of the nucleotides
        * @param deleted_atoms_ids A vector containing the indices of the atoms deleted because of the bonding between the nucleotides
        * @param chain_index The index of the strand in the system
        *
        * @sa setupChain
        * @sa generateConformerData
        * @sa HelicalParameters
        * @sa ConformationSearch
        */
        void setCoordsForChain(double *xyz, double *conf, PNAB::HelicalParameters &hp, std::vector<unsigned> &num_bu_atoms,
                               std::vector<unsigned> &bb_start_index, std::vector<double *> &base_coords_vec,
                               std::vector<unsigned> &deleted_atoms_ids, unsigned chain_index);
    };

}

#endif //PNAB_CHAIN_H
