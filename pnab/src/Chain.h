#ifndef PNAB_CHAIN_H
#define PNAB_CHAIN_H

#include <openbabel/forcefield.h>
#include <openbabel/rotor.h>
#include "Containers.h"
#define KJ_TO_KCAL 0.239006

namespace PNAB {

    /**
    * @brief A class for building nucleic acid strands and evaluating their energies
    *
    * The class creates the strands by connecting the nucleotides created in the PNAB::BaseUnit
    * class. It forms bonds between the nucleotides for the given sequence of the nucleobases,
    * and creates duplex or hexad systems if requested. The class also contains functions
    * for computing the energy terms for the system and checking whether the candidates are
    * accepted. Several functions in the class identifies the proper energy terms to be computed
    * for the system (e.g. bond and angle energies).
    *
    * @sa PNAB::ConformationSearch
    * @sa PNAB::HelicalParameters
    * @sa PNAB::RuntimeParameters
    * @sa PNAB::BaseUnit
    * @sa PNAB::Bases
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
        * @param strand A vector of the names of all bases in the strand, PNAB::RuntimeParameters::strand
        * @param ff_type The force field type, PNAB::RuntimeParameters::ff_type
        * @param range Monomer backbone index range
        * @param double_stranded Build double stranded system? PNAB::RuntimeParameters::is_double_stranded
        * @param hexad Build hexad system? PNAB::RuntimeParameters::is_hexad
        * @param strand_orientation The orientation of each strand in the hexad, PNAB::RuntimeParameters::strand_orientation
        *
        * @sa setupChain
        * @sa setupFFConstraints
        */
        Chain(PNAB::Bases bases, const PNAB::Backbone &backbone, std::vector<std::string> strand,
              std::string ff_type, std::array<unsigned, 2> &range, bool double_stranded = false, bool hexad = false,
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
        * @brief function to generate structure and energy data for nucleic acid conformers
        *
        * This func
        *
        * @sa Chain::setCoordsForChain
        * @sa Chain::fillConformerEnergyData
        * @sa ConformationSearch::RandomSearch
        * @sa ConformationSearch::MonteCarloSearch
        * @sa ConformationSearch::SystematicSearch
        * @sa ConformationSearch::GeneticAlgorithmSearch
        */ 
        PNAB::ConformerData generateConformerData(double *xyz, PNAB::HelicalParameters &hp, std::vector<double> energy_filter);
    
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
        std::vector<std::vector<std::vector<unsigned>>> v_fixed_bonds = std::vector<std::vector<std::vector<unsigned>>>(6);
        std::vector<std::vector<std::size_t>> v_base_connect_index = std::vector<std::vector<size_t>>(6);
        std::vector<std::vector<double*>> v_base_coords_vec_ = std::vector<std::vector<double*>>(6);
        unsigned chain_length_, n_chains_;
        bool isKCAL_, double_stranded_, hexad_, canonical_nucleobase_;
        std::vector<bool> strand_orientation_;
        OpenBabel::OBForceField *pFF_;
        std::array<unsigned, 2> monomer_bb_index_range_;
        std::vector<std::vector<unsigned>> v_bb_start_index_ = std::vector<std::vector<unsigned>>(6);
        std::string ff_type_;
        PNAB::Backbone backbone_;
        std::vector<std::vector<unsigned int>> all_torsions_;
        std::vector<bool> is_fixed_bond;
    
        OpenBabel::OBFFConstraints constraintsAng_, constraintsBond_, constraintsTor_, constraintsTot_;
    
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
        * @param deleted_atom_ids An empty vector to be populated with the indices of the atoms deleted because of the bonding between the nucleotides
        * @param num_base_unit_atoms An empty vector to be populated with the number of atoms in each BaseUnit
        * @param bb_start_index An empty vector to be populated with the starting index of the backbone atoms in the BaseUnit
        * @param base_coords_vec An empty vector to be populated with the coordinates of each one of the nucleotides
        * @param fixed_bonds_vec An empty vector to be populated with the indices of fixed bonds in each one of the nucleotides
        * @param chain_index The index of the chain
        *
        * @sa Chain::Chain
        */
        void setupChain(std::vector<PNAB::Base> &strand, OpenBabel::OBMol &chain, std::vector<unsigned> &new_bond_ids,
                            std::vector<unsigned> &deleted_atoms_ids, std::vector<unsigned> &num_base_unit_atoms,
                            std::vector<unsigned> &bb_start_index, std::vector<double *> &base_coords_vec, 
                            std::vector<std::vector<unsigned>> &fixed_bonds_vec, unsigned chain_index);
        /**
        * @brief Determines the terms that should be ignored during the computation of the bond, angle, and torsional energies
        *
        * For the bond energy, we consider only the one new bond that is formed between the adjacent nucleobases in the system.
        */
        void setupFFConstraints(OpenBabel::OBMol &chain, std::vector<unsigned> &new_bond_ids, std::vector<std::vector<unsigned>> &fixed_bonds_vec, unsigned offset = 0);
        void setCoordsForChain(double *xyz, double *conf, PNAB::HelicalParameters &hp, std::vector<unsigned> &num_bu_atoms,
                               std::vector<unsigned> &bb_start_index, std::vector<double *> &base_coords_vec,
                               std::vector<unsigned> &deleted_atoms_ids, unsigned chain_index);
    };
    
}

#endif //PNAB_CHAIN_H
