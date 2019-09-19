#ifndef PNAB_CONFORMATIONSEARCH_H
#define PNAB_CONFORMATIONSEARCH_H

#include "Chain.h"

#define BOLTZMANN 0.0019872041 // kcal/(mol.K)

namespace PNAB {
    /**
     * @brief A rotor search function used to find acceptable conformations of arbitrary backbone and helical
     * parameter combinations. The main class of the proto-Nucleic Acid Builder
     *
     * The general algorithm follows the following steps:
     * 1. Define the input parameters:
     *      - The three-dimensional structure of the backbone, the atoms that form the connection to the base,
     *        and the two atoms (head and tail) that link the backbone with two adjacent backbones.
     *      - The helical configuration of the nucleobases.
     *      - The runtime parameters:
     *          - The specific search algorithm (see below).
     *          - The distance and energy thresholds (see below).
     *          - The sequence of the nucleobases and whether the system is single-stranded, a duplex, or a hexad.
     *          .
     *      .
     * 2. Connect the backbone with the first nucleobase in the sequence.
     * 3. Fix all the bond distances, bond angles, and the rigid dihedral angles in the backbone.
     *    Fix all the atoms of the nucleobase.
     * 4. Change the dihdedal angles in the backbone using one of the search algorithm.
     * 5. Translate and rotate the tail atom of the backbone using the given helical parameters.
     * 6. Measure the distance between the head and the tail atoms of the backbone.
     * 7. If the distance is greater than the distance threshold, go to 4. Else, proceed.
     * 8. Arrange the nucleobases using the given helical parameters. Replicate the backbone structure
     *    across all the nucleotides.
     * 9. Evaluate various energy terms for the system.
     * 10. If the energy terms are less than the energy thresholds, save the structure. Else, reject the structure.
     * 11. Repeat steps 4-10 or quit.
     *
     * The various search algorithm differ in the implementation of step 4.
     * Four classes of search algorithms are implemented here, namely:
     * - Systematic Search
     * - Random Search
     * - Monte Carlo Search
     * - Genetic Algorithm Search
     *
     * Details on the specific implementation of each algorithm are discussed below.
     *
     * @sa SystematicSearch
     * @sa RandomSearch
     * @sa MonteCarloSearch
     * @sa GeneticAlgorithmSearch
     * @sa HelicalParameters
     * @sa RuntimeParameters
     * @sa Backbone
     * @sa BaseUnit
     * @sa Chain
     */
    class ConformationSearch {

    public:
        /** 
         * @brief Constructor for the conformation search class
         *
         * The constructor takes the input parameters and constructs a single nucleotide.
         * Then, it determines the rotatable dihedral angles in the backbone. 
         *
         * @param runtime_params Series of parameters that control how the search algorithm runs
         * @param backbone The backbone molecule over which the algorithm searches
         * @param helical_params The geometric parameters constraining the possible conformations of the backbone
         * @param bases List of defined bases that serves of as the library used for building the strand
         * @param prefix A string to prepend the name of the accepted backbone candidates
         * @param verbose Whether to print progress report to screen
         */
        ConformationSearch(PNAB::RuntimeParameters &runtime_params, PNAB::Backbone &backbone,
                           PNAB::HelicalParameters &helical_params, PNAB::Bases bases,
                           std::string prefix = "test", bool verbose = true);

        /**
         * @brief A function to call the appropriate search algorithm using the provided RuntimeParameters::search_algorithm
         *
         * @returns A string in CSV format containing the properties of the accepted candidates. The structure of the accepted
         *          candidates are saved as PDB files.
         *
         * @sa SystematicSearch
         * @sa RandomSearch
         * @sa MonteCarloSearch
         * @sa GeneticAlgorithmSearch
         */
        std::string run();

    private:
        bool verbose_; //!< @brief Whether to print progress report to screen
        PNAB::RuntimeParameters runtime_params_; //!< @brief The runtime parameters instance, RuntimeParameters
        std::array<unsigned, 2> backbone_range_; //!< @brief The Backbone index range for the first nucleotide
        PNAB::Backbone backbone_; //!< @brief The backbone molecule
        PNAB::HelicalParameters helical_params_; //!< @brief the helical parameters
        PNAB::Bases bases_; //!< @brief the list of the defined bases
        std::mt19937_64 rng_; //!< @brief A random number generator
        std::vector<PNAB::ConformerData> conf_data_vec_; //!< @brief A vector containing the data for the accepted candidates
        OpenBabel::matrix3x3 step_rot_, //!< @brief The step rotation matrix, HelicalParameters::getStepRotationOBMatrix
                             glbl_rot_; //!< @brief The global rotation matrix, HelicalParameters::getGlobalRotationOBMatrix
        OpenBabel::vector3 step_translate_, //!< @brief The step translation vector, HelicalParameters::getStepTranslationVec
                           glbl_translate_; //!< @brief The global translation vector, HelicalParameters::getGlobalTranslationVec
        OpenBabel::OBConversion conv_; //!< @brief An openbabel conversion object used to write accepted candidates
        OpenBabel::OBMol test_chain_; //!< @brief A molecule containing the structure of accepted candidates for saving
        unsigned monomer_num_coords_; //!< @brief 3*The number of atoms in the first nucleotide 
        std::string prefix_; //!< @brief A string to prepend the name of the accepted backbone candidates

        BaseUnit unit; //!< @brief The nucleotide unit for the first nucleotide in the strand. Used for searching for conformations
        OpenBabel::OBMol bu_a_mol; //!< @brief The molecule of the first nucleotide
        unsigned head, //!< @brief The first terminal atom in the backbone
                 tail; //!< @brief The second terminal atom in the backbone
        std::string output_string; //!< @brief A string in CSV format containing the properties of the accepted candidates
        double* coords; //!< @brief The coordinates of the first nucleotide 
        OpenBabel::OBRotorList rl; //!< @brief The list of the all rotatable dihedral angles in the backbone
        std::vector<OpenBabel::OBRotor*> rotor_vector; //!< @brief A vector of dihedral angles to be rotated in the search. Excludes fixed angles

        /**
         * @brief Given a step size, the algorithm exhaustively searches over all the rotatable dihedral angles in the backbone
         *
         * The systematic search algorithm works by exhaustively rotating all the dihedral angles between
         * 0 and 360&deg; using the given step size. This is the simplest search algorithm. However, it is
         * the most time-consuming as the time grows exponentially with the number of rotatable bonds. 
         * The number of steps required for the search is \f$(\frac{360.0}{\textrm{dihedral step}})^{\textrm{number of rotatable dihedrals}}\f$. 
         * This algorithm is deterministic and reproducible. It is guranteed to find an acceptable candidate, if any exists, to within
         * the given resolution. This algorithm can be used to validate the results of the other algorithms and to determine
         * whether the other algorithms found all the families of the acceptable candidates.
         * 
         * @sa RuntimeParameters::dihedral_step
         * @sa measureDistance
         * @sa Chain::generateConformerData
         */ 
        void SystematicSearch();

        /**
         *@brief The algorithm randomly changes all the dihedral angles in the backbone and evaluates whether they are acceptable
         *
         * At every iteration, this algorithm sets all the dihedral angles in the backbone to random values. This algorithm
         * is purely random and each iteration is completely independent from the previous iteration. 
         * This algorithm can sample the dihedral angle space much faster than the Systemtic Search algorithm. Therefore, it may
         * find an acceptable candidate faster. However, it is not deterministic.
         * As each iteration is completely independent, this algorithm does not become trapped in unfavorable configurations.
         * Nevertheless, the search does not improve with time. The algorithm finishes after the specified number of iterations. 
         *
         * The probability of choosing a random angle between 0 and 360&deg; can be uniform or can be weighted. See the
         * description in WeightedDistributions for an explanation on the weighting scheme for the dihedral angles.
         *
         * @param weighted Whether to use a weighted distribution for the dihedral angles
         *
         * @sa RuntimeParameters::num_steps
         * @sa RuntimeParameters::weighting_temperature
         * @sa WeightedDistributions
         * @sa measureDistance
         * @sa Chain::generateConformerData
         */ 
        void RandomSearch(bool weighted);

        /**
         * @brief The algorithm utilizes the Metropolis Monte Carlo scheme to improve the choice of the dihedral angle
         *
         * At every iteration, this algorithm randomly chooses two dihedrals and set them to random values.
         * If the new configuration decreases the distance between the head atom of the backbone and the tail atom
         * of the adjacent backbone, this step is accepted. If the new configuration does not decrease the distance,
         * the step is provisionally accpeted if a random number between 0 and 1 is less than 
         * \f$exp(\frac{E}{k_{B}T})\f$, where \f$E\f$ is a penalty term based on the distance, and \f$k_{B}T\f$ is the Boltzmann
         * constant multiplied by the temperature. If the step is provisionally accepted, the system is constructed
         * and the energy terms of the system are computed. If the energy terms are higher than the thresholds,
         * the step is finally rejected. The algorithm finishes after the specified number of iterations.
         *
         * The penalty term \f$E\f$ is given by \f$E=-k(\textrm{currrent distance} - \textrm{previous distance})^2\f$.
         * \f$k\f$ is a bond stretching constant fixed at 300 kcal/mol/Angstrom\f$^2\f$. The temperature controls how
         * agressively the algorithm should provisionally accept or reject steps. An infinite temperature will reduce
         * the algorithm to a random search with two dihedral angles changed at every step.
         *
         * The algorithm significantly accelerates the finding of dihedral angles that satisfy the distance threshold.
         * However, the distance threshold only indicates that the candidate backbone configuration can form a periodic
         * structure. It does not indicates whether the backbone configuration is low or high in energy. Therefore,
         * the ultimate criteria for accepting or rejecting a step, if the distance is not decreases in the step,
         * is whether it satisfies the energy thresholds.
         *
         * The probability of choosing a random angle between 0 and 360&deg; can be uniform or can be weighted. See the
         * description in WeightedDistributions for an explanation on the weighting scheme for the dihedral angles.
         * @param weighted Whether to use a weighted distribution for the dihedral angles
         *
         * @sa RuntimeParameters::num_steps
         * @sa RuntimeParameters::weighting_temperature
         * @sa RuntimeParameters::monte_carlo_temperature
         * @sa WeightedDistributions
         * @sa measureDistance
         * @sa Chain::generateConformerData
         */ 
        void MonteCarloSearch(bool weighted);

        /**
         * @brief This algorithm utilizes the genetic algorihtm procedure to improve the choice of the dihedral angle
         *
         * In the genetic algorithm search, a population of rotamers is initialized with random dihedral angles.
         * At each iteration (or generation), certain rotamers (parents) are chosen for reproduction based on their
         * fitness. The offsprings are generated by allowing the two parents to exchange (crossover) one dihedral angle.
         * Then, a random mutation (or a new value for the dihedral angle) may be introduced for one dihedral angle.
         * The new population consisting of the offsprings then follow the same procedure. The algorithm finishes
         * after the specified number of generations. The length of the search is proportional to the number of
         * generations multiplied by the size of the population.
         *
         * The fitness of the individuals is computed as the inverse of the distance between the head atom of the backbone
         * and the tail atom of the adjacent backbone. The individuals in the population are ranked by their fitness
         * score, and the probability of choosing individuals for mating is proportional to their ranking.
         * 
         * The exchange of the dihedral angles between the parents and the introduction of new dihedral angles
         * is determined by the crossover and mutation probabilities.
         *
         * The algorithm significantly accelerates the finding of dihedral angles that satisfy the distance threshold.
         * However, the distance threshold only indicates that the candidate backbone configuration can form a periodic
         * structure. It does not indicates whether the backbone configuration is low or high in energy. Therefore,
         * introducing mutations with higher probability is useful in preventing the formation of a homogeneous
         * population that has low distances but high energies.
         *
         * The probability of choosing a random angle between 0 and 360&deg; is uniform.
         *
         * @sa RuntimeParameters::num_steps
         * @sa RuntimeParameters::population_size
         * @sa RuntimeParameters::crossover_rate
         * @sa RuntimeParameters::mutation_rate
         * @sa measureDistance
         * @sa Chain::generateConformerData
         */ 
        void GeneticAlgorithmSearch();

        /**
         * @brief Produces weighted distributions for each rotatable dihedral angle in the backbone
         *
         * Instead of using a uniform distrbution for the dihedral angles between 0 and 360&deg;, this
         * function generates weighted distributions for each rotatable dihedral angle in the backbone.
         * The energy of each dihedral angle in isolation is computed every 5&deg;. Then, the
         * probability for a given angle is computed as the \f$\frac{exp(\frac{E_i}{k_{B}T})}{\sum_i exp(\frac{E_i}{k_{B}T})}\f$,
         * where \f$E_i\f$ is the dihedral energy at angle \f$i\f$. The probability in each interval is the
         * linear interpolation of the probabilities at the limits of the interval.
         * The temperature controls how agressively the algorithm should weight the dihedral angles.
         * An infinite temperature will reduce the probabilities to uniformly random distributions.
         *
         * This weighting scheme is useful in increasing the sampling of the dihedral angles where
         * the energy of the dihedral angle is low. However, the energy of the dihedral angle in isolation
         * does not necessarily mean the nucleotide energy is going to be low. This weighting scheme may improve
         * the search if the acceptable candidates adopt structures where the dihedral angles are not strained. However,
         * it may worsen the search if the unstrained dihedral angles do not produce low-energy candidates.
         *
         * @returns A vector of the weighted distributions for each rotatable dihedral angle
         *
         * @sa RuntimeParameters::weighting_temperature
         * @sa RandomSearch
         * @sa MonteCarloSearch
         */ 
        std::vector <std::piecewise_linear_distribution<double>> WeightedDistributions();

        /**
        * @brief Compute the distance between the head backbone atom and the tail backbone atom of the next nucleotide
        *
        * This function applies the global translation and rotation to both the head and the tail atoms of the first nucleotide.
        * Then, it applies the step translation and rotation to the tail atom to get its coordinates in the second nucleotide.
        * Then, it computes the distance between the two atoms. This function is called by all the search algorithms.
        *
        * @param coords The coordinates of the first nucleotide
        * @param head The index of the first terminal atom in the backbone of the first nucleotide
        * @param tail The index of the second terminal atom the backbone of the first nucleotide
        *
        * @returns The distance between the two atoms in Angstroms
        *
        * @sa RuntimeParameters::max_distance
        * @sa ConformationSearch
        */ 
        double measureDistance(double *coords, unsigned head, unsigned tail);

        /**
        * @brief A function to report the data on the accepted candidates 
        *
        * This function saves the structure of each accepted candidate in PDB format. It also reports
        * the properties of the candidates, ordered by their total energies, in a string with the CSV format.
        * This is called by all search algorithms when a candidate is accepted. The ConformationSearch::output_string
        * variable is updated.
        *
        * @param conf_data The properties of the accepted candiate
        */ 
        void reportData(PNAB::ConformerData conf_data);

        /**
        * @brief Calculate the RMSD between a given accpeted candidate the accepted candidate with the lowest energy
        *
        * @param ref The coordinates of the lowest energy candidate
        * @param conf The coordinates of the accepted candidate
        * @param size The number of atoms*3 in the molecule
        *
        * @returns The RMSD value
        *
        * @sa reportData
        */ 
        double calcRMSD(double *ref, double *conf, unsigned long size) {
            double rmsd = 0;
            for (int i = 0; i < size; i++) {
                rmsd += pow(*(ref + i) - *(conf + i),2);
            }
            return sqrt(rmsd/(size/3));
        }

        /**
        * @brief Prints the percentage of search completed and the best accepted candidate
        *
        * Prints to the standard output the percentage of the search completed, the name of the PDB file
        * containing the best accepted candidte, its distance and energy, and the number of accepted candidates.
        *
        * @param search_index The step number in the search
        * @param search_size The total number of steps
        */ 
        void printProgress(std::size_t search_index, std::size_t search_size);
    };
}

#endif //PNAB_MONTECARLOROTORSEARCH_H
