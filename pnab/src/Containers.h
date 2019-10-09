/**@file
 * @brief A file for declaring various classes for defining options
 */

#ifndef PNAB_CONTAINERS_H
#define PNAB_CONTAINERS_H

#include <string>
#include <random>
#include <array>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

//! The PNAB name space contains all the C++ classes and functions for the proto-Nucleic Acid Builder.
namespace PNAB {
    /**
    * @brief A class for holding necessary and optional runtime parameters for conformational searches
    *
    * The runtime parameters are used for building the strands and during the conformational search.
    *
    * @sa Chain
    * @sa ConformationSearch
    */
    class RuntimeParameters {

    public:
        /**
        * @brief Empty constructor.
        *
        * This empty constructor can be used. After that, values for the member variables should be specified.
        */
        RuntimeParameters() : energy_filter{}, max_distance(), ff_type(),
                              num_steps(0), seed(0), weighting_temperature(298.0), monte_carlo_temperature(298.0),
                              mutation_rate(0.75), crossover_rate(0.75), population_size(1000), strand{}, is_hexad(false),
                              build_strand(std::vector<bool> {true, false, false, false, false, false}),
                              strand_orientation(std::vector<bool> {true, true, true, true, true, true}){};

        // Thresholds
        std::vector<double> energy_filter;      /*!< @brief [max bond E, max angle E, max torsion E, max VDW E, max total E]
                                                *
                                                * - Maximum accepted energy for newly formed bonds in the backbone (kcal/mol/bond)
                                                * - Maximum accepted energy for newly formed angles in the backbone (kcal/mol/angle)
                                                * - Maximum accepted torsional energy for rotatable bonds (kcal/mol/nucleotide)
                                                * - Maximum accepted van der Waals energy (kcal/mol/nucleotide)
                                                * - Maximum accepted total energy (kcal/mol/nucleotide)
                                                *
                                                * @sa Chain::generateConformerData
                                                * @sa Chain::fillConformerEnergyData
                                                */

        double max_distance;                    /*!< @brief Maximum accepted distance (Angstrom) between head and tail of successive nucleotides
                                                *
                                                * This distance is used to quickly screen rotamers. The algorithm searches for backbone
                                                * candidates that are periodic in terms of the helical structure. Thus, the terminal atoms
                                                * in adjacent candidates must be within a small distance to allow for a bond to form.
                                                * This distance is recommended to be less than 0.1 Angstroms. The extra terminal atom
                                                * in the adjacent nucleotide is then removed.
                                                *
                                                * @sa ConformationSearch::measureDistance
                                                */

        // Force Field Parameter
        std::string ff_type;                    /*!< @brief The type of the forcefield such as "GAFF" or "MMFF94"; available through Openbabel
                                                *
                                                * @sa Chain::Chain
                                                */

        // Algorithm parameters
        std::string search_algorithm;           /*!< @brief The search algorithm
                                                *
                                                * There are six search algorithms:
                                                * - Systematic search
                                                * - Monte Carlo search
                                                * - Weighted Monte Carlo search
                                                * - Random search
                                                * - Weighted random search
                                                * - Genetic algorithm search
                                                *
                                                *   @sa ConformationSearch
                                                */

        std::size_t num_steps;                  /*!< @brief The number of points sampled in Monte Carlo and random searches
                                                * and the number of generations in the genetic algorithm search
                                                *
                                                * @sa ConformationSearch::MonteCarloSearch
                                                * @sa ConformationSearch::RandomSearch
                                                * @sa ConformationSearch::GeneticAlgorithmSearch
                                                */
        unsigned int seed;                      /*< @brief The seed for the random number generator. Use the same value to get reproducible
                                                * results.
                                                *
                                                *  @sa ConformationSearch::MonteCarloSearch
                                                *  @sa ConformationSearch::RandomSearch
                                                *  @sa ConformationSearch::GeneticAlgorithmSearch
                                                */

        double dihedral_step;                   /*!< @brief The dihedral step size for systematic search (degrees)
                                                *
                                                * The number of steps will be \f$(\frac{360.0}{\textrm{dihedral step}})^{\textrm{number of rotatable dihedrals}}\f$.
                                                *
                                                * @sa ConformationSearch::SystematicSearch
                                                */

        double weighting_temperature;           /*!< @brief The temperature used to compute the weighted probability for weighted
                                                * Monte Carlo and weighted random searches
                                                *
                                                * @sa ConformationSearch::WeightedDistributions
                                                */

        double monte_carlo_temperature;         /*!< @brief The temperature used in the Monte Carlo acceptance and rejection procedure
                                                *
                                                * @sa ConformationSearch::MonteCarloSearch
                                                */

        double mutation_rate;                   /*!< @brief The mutation rate in the genetic algorithm search
                                                *
                                                * @sa ConformationSearch::GeneticAlgorithmSearch
                                                */

        double crossover_rate;                  /*!< @brief The crossover rate in the genetic algorithm search
                                                *
                                                * @sa ConformationSearch::GeneticAlgorithmSearch
                                                */

        int population_size;                    /*!< @brief The population size in the genetic algorithm search
                                                *
                                                * @sa ConformationSearch::GeneticAlgorithmSearch
                                                */

        //Strand parameters
        std::vector<std::string> strand;        //!< @brief The names of each base used in the strand
        std::vector<bool> build_strand;         /*!< @brief Defines whether to build a given strand
                                                *
                                                * @sa Chain::Chain
                                                */
        std::vector<bool> strand_orientation;   /*!< @brief Defines strand orientation for each strand in the hexad
                                                *
                                                * @sa Chain::setCoordsForChain
                                                * @sa Chain::setCoordsForChain
                                                */
        bool is_hexad;                          //!< @brief Defines whether the 60 degrees rotation for hexads is performed
    };

    /**
    * @brief A class for holding values for all helical parameters
    *
    * The helical parameters are used for generating the geometries of the nucleobases in the strands.
    * This class holds the value for six helical parameters and functions to generate the geometries.
    *
    * @sa Chain::setCoordsForChain
    * @sa ConformationSearch::measureDistance
    */
    class HelicalParameters {

    public:
        /**
        * @brief Empty constructor
        *
        * This empty constructor can be used. After that, values for the member variables should be specified.
        */
        HelicalParameters() : h_twist{0}, h_rise{0}, inclination{0}, tip{0}, x_displacement{0}, y_displacement{0} {};

        //Helical parameters
        double     inclination,                         //!< @brief Inclination
                   tip,                                 //!< @brief Tip
                   h_twist,                             //!< @brief Helical twist
                   x_displacement,                      //!< @brief X-Displacement
                   y_displacement,                      //!< @brief Y-Displacement
                   h_rise;                              //!< @brief Helical rise

        /**
        * @brief Get the global rotation matrix in the OpenBabel::matrix3x3 format
        *
        * @returns The global rotation matrix
        *
        * @sa getGlobalRotationMatrix
        */
        OpenBabel::matrix3x3 getGlobalRotationOBMatrix() {
            auto arr = getGlobalRotationMatrix();
            return OpenBabel::matrix3x3(OpenBabel::vector3(arr[0], arr[1], arr[2]),
                                        OpenBabel::vector3(arr[3], arr[4], arr[5]),
                                        OpenBabel::vector3(arr[6], arr[7], arr[8]));
        }

        /**
        * @brief Get the global translation vector using the HelicalParameters::x_displacement and HelicalParameters::y_displacement
        *
        * @returns The translation vector
        */
        OpenBabel::vector3 getGlobalTranslationVec() {
            return OpenBabel::vector3(x_displacement, y_displacement, 0);
        }

        /**
        * @brief Get the step rotation matrix in the OpenBabel::matrix3x3 format
        *
        * @param n The sequence of the nucleobase in the strand
        *
        * @returns The step rotation matrix
        *
        * @sa getStepRotationMatrix
        */
        OpenBabel::matrix3x3 getStepRotationOBMatrix(unsigned n = 0) {
            auto arr = getStepRotationMatrix(n);
            return OpenBabel::matrix3x3(OpenBabel::vector3(arr[0], arr[1], arr[2]),
                                        OpenBabel::vector3(arr[3], arr[4], arr[5]),
                                        OpenBabel::vector3(arr[6], arr[7], arr[8]));
        }

        /**
        * @brief Get the step translation vector using the HelicalParameters::h_rise
        *
        * @param n The sequence of the nucleobase in the strand
        *
        * @returns The step translation vector
        */
        OpenBabel::vector3 getStepTranslationVec(unsigned n = 0) {
            return OpenBabel::vector3(0, 0, n * h_rise);
        }

    private:
        /**
        * @brief Get the global rotation matrix using the HelicalParameters::tip and HelicalParameters::inclination
        *
        * @returns The global rotation matrix
        *
        * @sa rodrigues_formula
        */

        // Lu, X. J., El Hassan, M. A., & Hunter, C. A. (1997). Structure and conformation of helical nucleic acids:
        // rebuilding program (SCHNArP). Journal of molecular biology, 273(3), 681-691.
        std::array<double, 9> getGlobalRotationMatrix() {
            double eta = inclination * DEG_TO_RAD, theta = tip * DEG_TO_RAD;
            double Lambda = sqrt(eta*eta + theta*theta);
            
            std::array<double, 3> axis;
            if (Lambda != 0)
                axis = {eta/Lambda, theta/Lambda, 0};
            else
               axis = {0, 1, 0};

            return rodrigues_formula(axis, Lambda);
        }

        /**
        * @brief Get the step rotation matrix using the HelicalParameters::h_twist
        *
        * @param n The sequence of the nucleobase in the strand
        *
        * @returns The step rotation matrix
        *
        * @sa matrix_mult
        */
        std::array<double, 9> getStepRotationMatrix(unsigned n = 0) {
            double Omega = h_twist * DEG_TO_RAD;
            std::array<double, 9> m_mat {cos(Omega), -sin(Omega), 0, sin(Omega), cos(Omega), 0, 0, 0, 1};
            std::array<double, 9> r_mat = {1, 0, 0, 0, 1, 0, 0, 0, 1};
            for (int i = 0; i < n; ++i)
                r_mat = matrix_mult(m_mat, r_mat);

            return r_mat;
        };


        /**
        * @brief Rodrigues rotation formula for rotating a vector in space
        *
        * Outputs a 3x3 matrix in the form of a one dimensional array to be used
        *
        * @param axis A unit vector defining the axis about which to rotate by angle theta
        * @param theta The angle at which to rotate about vector given by axis
        *
        * @return The new rotation matrix
        */
        std::array<double, 9> rodrigues_formula(std::array<double, 3> axis, double theta) {
            std::array<double, 9> m{};
            m[0] = cos(theta) + axis[0]*axis[0]*(1 - cos(theta));
            m[1] = axis[0]*axis[1]*(1 - cos(theta)) - axis[2]*sin(theta);
            m[2] = axis[1]*sin(theta) + axis[0]*axis[2]*(1 - cos(theta));
            m[3] = axis[2]*sin(theta) + axis[0]*axis[1]*(1 - cos(theta));
            m[4] = cos(theta) + axis[1]*axis[1]*(1 - cos(theta));
            m[5] = -axis[0]*sin(theta) + axis[1]*axis[2]*(1 - cos(theta));
            m[6] = -axis[1]*sin(theta) + axis[0]*axis[2]*(1 - cos(theta));
            m[7] = axis[0]*sin(theta) + axis[1]*axis[2]*(1 - cos(theta));
            m[8] = cos(theta) + axis[2]*axis[2]*(1 - cos(theta));

            return m;
        };

        /**
        * @brief Matrix multiplication function
        *
        * @param m1 Matrix 1
        * @param m2 Matrix 2
        *
        * @returns The matrix product
        */
        std::array<double, 9> matrix_mult(std::array<double, 9> m1, std::array<double, 9> m2) {
            std::array<double, 9> out{};

            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    for (int k = 0; k < 3; ++k)
                        out[3*i + j] += m1[3*i + k]*m2[j + 3*k];

            return out;
        };
    };

    /**
    * @brief Class for holding backbone information
    *
    * The backbone here refers to the backbone in a single nucleotide. This class holds information
    * on the molecular structure of the backbone and the bonds that the backbone form with the
    * nucleobases and the adjacent backbones. This class also has functions to manipulate the backbone.
    *
    * @sa BaseUnit
    * @sa ConformationSearch
    */
    class Backbone {
    public:

        /**
        * @brief Empty constructor
        *
        * This empty constructor can be used. After that, values for the member variables should be specified.
        */
        Backbone() : backbone{}, file_path{}, interconnects{}, linker{}, vector_atom_deleted{}, fixed_bonds{} {}

        /**
        * @brief Constructor for the backbone unit
        *
        * @param file_path The path to the file containing the molecule.
        * @param interconnects The atom indices that define the periodic conditions between backbones { head, tail }.
        * @param linker The atom indices used to align and connect backbone to base in the nucleotide.
        * @param fixed_bonds A vector containing pairs of indices defining fixed rotatable bonds during dihedral search.
        */
        Backbone(std::string file_path, std::array<unsigned, 2> interconnects, std::array<unsigned,2> linker, std::vector<std::vector<unsigned>> fixed_bonds = {});

        /**
        * @brief Gives the pointer to an atom that is the head from Backbone::interconnects{head, tail}
        * @return The atom pointer from the backbone OBMol object
        */
        OpenBabel::OBAtom* getHead() {
            return backbone.GetAtom(interconnects[0]);
        }

        /**
        * @brief Gives the pointer to an atom that is the tail from Backbone::interconnects{head, tail}
        * @return Pointer to the atom for the tail
        */
        OpenBabel::OBAtom* getTail() {
            return backbone.GetAtom(interconnects[1]);
        }

        /**
        * @brief Get the first Backbone::linker atom pointer
        * @return Pointer to the atom that is the one linking to the Base
        */
        OpenBabel::OBAtom* getLinker() {
            return backbone.GetAtom(linker[0]);
        }

        /**
        * @brief Get the second Backbone::linker atom pointer (which is probably a hydrogen)
        * @return Pointer to the atom that defines the vector from the backbone to the base
        */
        OpenBabel::OBAtom* getVector() {
            if (!vector_atom_deleted)
                return backbone.GetAtom(linker[1]);
            else {
                std::cerr << "Called getVector() for backbone with no vector atom." << std::endl;
                return nullptr;
            }
        }

        /**
        * @brief Centers the molecule. Basically just an alias of the Center() function from OpenBabel.
        */
        void center() {
            backbone.Center();
        }

        /**
        * @brief Rotates the molecule by a matrix. Basically just an alias of the Rotate() function from OpenBabel.
        * @param rot The matrix by which to rotate the molecule
        */
        void rotate(double* rot) {
            backbone.Rotate(rot);
        }

        /**
        * @brief Translates the molecule by a vector. Basically just an alias of the Translate() function from OpenBabel.
        * @param vec The vector by which to translate the molecule
        */
        void translate(OpenBabel::vector3 vec) {
            backbone.Translate(vec);
        }

        /**
        * @brief Gives a copy of the molecule in the backbone, Backbone::backbone
        * @return A copy of the backbone molecule
        */
        OpenBabel::OBMol getMolecule() {
            return backbone;
        }

        /**
        * @brief Deletes the atom from getVector() safely. If the atom is already deleted, nothing happens.
        */
        void deleteVectorAtom();

        std::array<unsigned , 2> interconnects,          //!< @brief The atom indices that define the periodic conditions between backbones { head, tail }
                                 linker;                 //!< @brief The atom indices used to align and connect backbone to base in the nucleotide
        std::vector<std::vector<unsigned>> fixed_bonds;  //!< @brief A vector containing pairs of indices defining fixed rotatable bonds during dihedral search
        OpenBabel::OBMol backbone;                       //!< @brief The molecule for the backbone
        std::string file_path;                           //!< @brief The path to the file containing the molecule

    private:
       /**
       * @brief Does some basic sanity checks (such as whether or not the indices of the atom are within the range
       * of the molecule)
       */
        void validate();

        bool vector_atom_deleted;                        //!< @brief Whether or not the atom from getVector() has been deleted
    };

    /**
    * @brief Class to fully define bases (i.e. Adenine, Cytosine)
    *
    * This class holds information on the molecular structure of one nucleobase
    * and the bond that it should form with the backbone. This class also has functions
    * to manipulate the nucleobase and gets information about it.
    *
    * @sa Bases
    * @sa BaseUnit
    */
    class Base {

    public:

        /**
        * @brief Empty constructor
        *
        * This empty constructor can be used. After that, values for the member variables should be specified.
        */
        Base() : name{}, code{}, linker{}, base{}, vector_atom_deleted{}, pair_name{} {};

        /**
        * @brief Create Base from basic set of parameters
        * @param name name of the base
        * @param code three-letter code
        * @param file_path path to the backbone file
        * @param linker indices for atoms forming the vector connecting to the backbone
        * @param pair_name Name of the pairing base
        */

        Base(std::string name, std::string code, std::string file_path, std::array<std::size_t, 2> linker,
             std::string pair_name = "");

        /**
        * @brief Gives the atom of the base that connects directly to the backbone, Base::linker[0]
        * @return A pointer to the atom that connects to the backbone
        */
        OpenBabel::OBAtom* getLinker() {
            return base.GetAtom(static_cast<unsigned>(linker[0]));
        }

        /**
        * @brief Gives the (most likely hydrogen) atom of the base connected to the atom from getLinker() which defines how
        * the base connects, Base::linker[1]
        * @return A pointer to the atom forming the vector connecting to the backbone
        */
        OpenBabel::OBAtom* getVector() {
            if (!vector_atom_deleted)
                return base.GetAtom(static_cast<unsigned>(linker[1]));
            else
                return nullptr;
        }

        /**
        * @brief Returns a copy of the base molecule, Base::base
        * @return A copy of the base molecule
        */
        OpenBabel::OBMol getMolecule() {
            return base;
        }

        /**
        * @brief Gives the three-letter code of the base, Base::code
        * @return The code of the base
        */
        std::string getCode() {
            return code;
        }

        /**
        * @brief Gives the full name of the base, Base::name
        * @return The full name of the base
        */
        std::string getName() {
            return name;
        }

        /**
        * @brief Deletes the atom from getVector() safely. If the atom is already deleted, nothing happens.
        */
        void deleteVectorAtom() {
            if (!vector_atom_deleted) {
                size_t id = getLinker()->GetId();
                base.DeleteAtom(getVector());
                linker = {base.GetAtomById(id)->GetIdx(), 0};
                vector_atom_deleted = true;
            }
        }

        /**
        * @brief Get the name of the pair base, Base::pair_name
        * @return Pair name
        */
        std::string getBasePairName() {
            return pair_name;
        }

        std::string name,                               //!< @brief Full name of base (i.e. "Adenine" or just "A")
                    code,                               //!< @brief Three character code to define base ("Adenine": "ADE")
                    pair_name,                          //!< @brief Name of the pair base
                    file_path;                          //!< @brief Path to a file containing the base
        OpenBabel::OBMol base;                          //!< @brief The OBMol defining the base
        std::array<std::size_t, 2 > linker;             //!< @brief Holds indices for atoms forming a vector to connect to backbone {linker, hydrogen}

    private:
        /**
        * @brief Does some basic sanity checks (such as whether or not the indices of the atom are within the range
        * of the molecule).
        */
        void validate();

        bool vector_atom_deleted;                       //!< @brief Whether or not the getVector() atom was deleted
    };

    /**
    * @brief A class that contains a vector of all the defined bases and a funtion to
    * return a base and the complimentary base for all the bases defined.
    *
    * @sa Base
    */
    class Bases {
    public:
        /**
        * @brief Basic constructor for the Bases.
        *
        * The given bases are processed to create a vector of the defined bases and a vector of the complimentary bases.
        * This calls Base to make sure the given bases have Openbabel molecules
        *
        * @param input_bases The vector containing the information about the bases needed
        */
        Bases(std::vector<Base> input_bases);

        /**
        * @brief Empty constructor.
        */
        Bases() {};

        /**
        * @brief Returns the Base instance given the name of the base.
        *
        * @param name name of the base
        *
        * @returns the base
        */
        PNAB::Base getBaseFromName(std::string name) {
            for (auto v : bases) {
                if (v.getName().find(name) != std::string::npos)
                    return v;
            }
            std::cerr << "Base \"" << name << "\" does not exists in list of bases. Please check input file."
                      << std::endl;
            throw 1;
        }

        /**
        * @brief Returns the vector of the instances of Base given the names of the bases in the strand.
        *
        * @param strand vector of base names in the strand
        *
        * @returns vector of bases in the strand
        *
        * @sa Chain::Chain
        */
        std::vector<Base> getBasesFromStrand(std::vector<std::string> strand);

        /**
        * @brief Returns the complimentary vector of the instances of Base given the names of the bases in the strand.
        *
        * @param strand vector of base names in the strand
        *
        * @returns vector of bases in the complimentary strand
        *
        * @sa Chain::Chain
        */
        std::vector<Base> getComplimentBasesFromStrand(std::vector<std::string> strand);

    private:
        std::vector<Base> bases;                          //!< @brief The vector of bases
        bool all_bases_pair;                              //!< @brief Whether all the bases in the strand have complimentary bases
        std::map<std::string, PNAB::Base> name_base_map;  //!< @brief A map of the names of the bases and the complimentary bases
    };

    /**
    * @brief Class to hold bases with backbones attached (nucleotides), along with associated necessary information
    *
    * @sa Base
    * @sa Backbone
    * @sa Chain::setupChain
    * @sa ConformationSearch
    */
    class BaseUnit {
    public:
        /**
        * Constructor for the base unit
        *
        * @param b Base instance
        * @param backbone Backbone instance
        */
        BaseUnit(Base b, Backbone backbone);

        /**
        * @brief Empty constructor.
        */
        BaseUnit() {};

        /**
        * @brief Returns the nucleotide molecule, BaseUnit::unit
        *
        * @return The nucleotide molecule
        */
        const OpenBabel::OBMol getMol() {
            return unit;
        }

        /**
        * @brief Returns the indices for the begining and end of the nucleobase atom indices, BaseUnit::base_index_range
        *
        * @return The indices for the limits of the nucleobases
        */
        const std::array< std::size_t, 2 > getBaseIndexRange() {
            return base_index_range;
        };

        /**
        * @brief Returns the indices for the begining and end of the backbone atom indices, BaseUnit::backbone_index_range
        *
        * @return The indices for the limits of the backbone
        */
        const std::array< std::size_t, 2 > getBackboneIndexRange() {
            return backbone_index_range;
        };

        /**
        * @brief Returns the indices for atoms where the backbone connects, BaseUnit::backbone_interconnects
        *
        * @return The indices for atoms where the backbone connects
        */
        const std::array< std::size_t, 2 > getBackboneLinkers() {
            return backbone_interconnects;
        };

        /**
        * @brief Returns the index of the atom where the nucleobase connects to the backbone, BaseUnit::base_connect_index
        *
        * @return The nucleobase index connecting to the backbone
        */
        std::size_t getBaseConnectIndex() {
            return base_connect_index;
        }

        /**
        * @brief Returns a vector of the pair of indices for fixed rotatable dihedrals in the backbone, BaseUnit::fixed_bonds
        *
        * @return The indices of the fixed rotatable dihedrals
        */
        std::vector<std::vector<unsigned>> getFixedBonds() {
            return fixed_bonds;
        }

    private:
        OpenBabel::OBMol unit;                                                       //!< @brief Holds molecule containing base with backbone attached
        std::array< std::size_t, 2 > base_index_range,                               //!< @brief Range of indices of the unit that are a part of the base, [start, stop]
                                     backbone_index_range;                           //!< @brief Range of indices of the unit that are a part of the backbone, [start, stop]
        std::size_t base_connect_index;                                              //!< @brief Atom index where backbone connects to base (the base atom)
        std::array< std::size_t, 2 > backbone_interconnects;                         //!< @brief Atom indices defining where backbone connects
        std::vector<std::vector<unsigned>> fixed_bonds;                              //!< @brief Indices of fixed bonds in dihedral search
    };

    /**
     * @brief Class to contain important information for an individual conformer
     *
     * It includes detailed information about the energy components important
     * for distinguishing between different conformers. It also includes the value
     * of the RMSD relative to the best candidate. If the conformer satisfies
     * the distance and energy thresholds, then it is saved.
     *
     * @sa Chain::generateConformerData
     * @sa ConformationSearch::reportData
     * @sa RuntimeParameters::energy_filter
     * @sa RuntimeParameters::max_distance
     */
    struct ConformerData {
        double *coords,                             //!< @brief Pointer to array containing coordinates of all atoms in molecule chain
               *monomer_coord,                      //!< @brief Pointer to array containing coordinates of a single monomer
                distance,                           //!< @brief distance between interconnects in Backbone for adjacent BaseUnit
                bondE,                              //!< @brief Energy of newly formed bonds in the backbone divided by the length of the strand -1
                angleE,                             //!< @brief Energy of newly formed angles in the backbone divided by the length of the strand -1
                torsionE,                           //!< @brief Energy of all rotatable torsions divided by the length of the strand
                VDWE,                               //!< @brief Total van Der Wals Energy divided by the length of the strand
                total_energy,                       //!< @brief Total energy of the conformation divided by divided by the length of the strand
                rmsd;                               //!< @brief Root-mean square distance relative to lowest energy conformer
        std::size_t index;                          //!< @brief The index of the conformer
        bool chain_coords_present;                  //!< @brief Have the chain coordinates in coord been deleted?
        bool accepted;                              //!< @brief Is the energy of the conformer less than the thresholds

        /**
         * @brief Used for simple sorting based on total energy of the conformer
         * @param cd The ConformerData element to compare current element to
         * @return True if the other ConformerData has greater total energy, false otherwise
         */
        bool operator<(const ConformerData &cd) const {
            return (total_energy < cd.total_energy);
        }
    };

}

#endif //PNAB_CONTAINERS_H
