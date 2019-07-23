//
// Created by jbarnett8 on 8/22/17.
//

#ifndef PNAB_CONTAINERS_H
#define PNAB_CONTAINERS_H

#include <string>
#include <random>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

namespace PNAB {
    /**
     * \brief A class for holding necessary parameters for conformational searches
     */
    class RuntimeParameters {

    public:
        /**
         * \brief Empty constructor. Should not generally be used.
         */
        RuntimeParameters() : energy_filter{}, max_distance(), type(),
                              num_steps(0), strand{}, is_double_stranded(false),
                              is_hexad(false), is_parallel(true){};

        // Energy parameters
        std::vector<double> energy_filter;      //!< \brief {max bond E, max angle E, max VDW E, max total E}
        double max_distance;                    //!< \brief The maximum distance between head and tail of successive UnitChains that is accepted

        // Force Field Parameters
        std::string type;                       //!< \brief The type of the ForceField such as "GAFF" or "MMFF94"

        std::size_t num_steps;                  //!< \brief Determines how many points are sampled in the Monte Carlo searches

        //Strand parameters
        std::vector<std::string> strand;        //!< \brief Holds the names of each base used in the strand
        bool is_double_stranded;                //!< \brief Defines if the program prints and tests double stranded molecules
        bool is_hexad;                          //!< \brief Defines if the program prints and tests hexad molecules
        bool is_parallel;                       //!< \brief Defines if the hexad molecule is parallel or anti-parallel

    };

    /**
     * \brief Holds values for all helical parameters
     */
    class HelicalParameters {

    public:
        HelicalParameters() : tilt{0}, roll{0}, twist{0}, shift{0}, slide{0}, rise{0}, buckle{0}, propeller{0},
                              opening{0}, shear{0}, stretch{0}, stagger{0}, inclination{0}, tip{0}, x_displacement{0},
                              y_displacement{0} {};

        std::array<double, 9> getGlobalRotationMatrix() {
            double eta = inclination * DEG_TO_RAD, theta = tip * DEG_TO_RAD;
            double phi_pp = atan(std::isnan(eta / theta) ? 0 : eta / theta), Lambda = sqrt(eta*eta + theta*theta);
            std::array<double, 3> axis{sin(phi_pp), cos(phi_pp), 0};

            return rodrigues_formula(axis, Lambda);
        }

        OpenBabel::matrix3x3 getGlobalRotationOBMatrix() {
            auto arr = getGlobalRotationMatrix();
            return OpenBabel::matrix3x3(OpenBabel::vector3(arr[0], arr[1], arr[2]),
                                        OpenBabel::vector3(arr[3], arr[4], arr[5]),
                                        OpenBabel::vector3(arr[6], arr[7], arr[8]));
        }

        OpenBabel::vector3 getGlobalTranslationVec() {
            return OpenBabel::vector3(x_displacement, y_displacement);
        }

        //std::array<double, 9> getBasePairRotationMatrix() {
        //};

        OpenBabel::vector3 getBasePairTranslationVec() {
            return OpenBabel::vector3(shear, stretch, stagger);
        }

        std::array<double, 9> getStepRotationMatrix(unsigned n = 0) {
            double tau = tilt * DEG_TO_RAD, rho = roll * DEG_TO_RAD, Omega = twist * DEG_TO_RAD;
            double phi = atan(std::isnan(tau / rho) ? 0 : tau / rho), Gamma = sqrt(tau*tau + rho*rho);
            std::array<double, 3> axis{sin(phi), cos(phi), 0};

            auto m = rodrigues_formula(axis, Gamma);

            std::array<double, 9> z{cos(Omega), -sin(Omega), 0, sin(Omega), cos(Omega), 0, 0, 0, 1};

            auto m_mat = matrix_mult(z, m);
            std::array<double, 9> r_mat = {1, 0, 0, 0, 1, 0, 0, 0, 1};
            for (int i = 0; i < n; ++i)
                r_mat = matrix_mult(m_mat, r_mat);
            return r_mat;
        };

        OpenBabel::matrix3x3 getStepRotationOBMatrix(unsigned n = 0) {
            auto arr = getStepRotationMatrix(n);
            return OpenBabel::matrix3x3(OpenBabel::vector3(arr[0], arr[1], arr[2]),
                                        OpenBabel::vector3(arr[3], arr[4], arr[5]),
                                        OpenBabel::vector3(arr[6], arr[7], arr[8]));
        }

        OpenBabel::vector3 getStepTranslationVec(unsigned n = 0) {
            //n++;
            return OpenBabel::vector3(n * shift, n * slide, n * rise);
        }

        double getTwist() {
            return twist;
        }

        //Step parameters (describe rotations and translations between successive base-pairs)
        double     tilt,                                //!< \brief Symbol: tau; step parameter; rotation in x axis
                   roll,                                //!< \brief Symbol: rho; step parameter; rotation in y axis
                   twist,                               //!< \brief Symbol: Omega; step parameter; rotation in z axis
                   shift,                               //!< \brief Symbol: Dx; step parameter; translation in x
                   slide,                               //!< \brief Symbol: Dy; step parameter; translation in y
                   rise;                                //!< \brief Symbol: Dz; step parameter; translation in z

        //Base-pair parameters (describe rotations and translations between two bases
        double     buckle,                              //!< \brief Symbol: kappa; base-pair parameter; rotation in x axis
                   propeller,                           //!< \brief Symbol: omega; base-pair parameter; rotation in y axis
                   opening,                             //!< \brief Symbol: sigma; base-pair parameter; rotation in z axis
                   shear,                               //!< \brief Symbol: Sx; base-pair parameter; translation in x axis
                   stretch,                             //!< \brief Symbol: Sy; base-pair parameter; translation in y axis
                   stagger;                             //!< \brief Symbol: Sz; base-pair parameter; translation in z axis

        //Global parameters (describe rotations and translations of every base-pair)
        double     inclination,                         //!< \brief Symbol: eta; global parameter; rotation in x axis
                   tip,                                 //!< \brief Symbol: theta; global parameter; rotation in y axis
                   // twist_g                           ignored, only shown for completeness
                   x_displacement,                      //!< \brief Symbol: eta; global parameter; translation in x axis
                   y_displacement;                      //!< \brief Symbol: eta; global parameter; translation in y axis
                   // rise_g                            ignored, only shown for completeness

        // Random number generator
        std::mt19937_64 rng;

        /**
         * \brief Outputs a 3x3 matrix in the form of a one dimensional array to be used
         * @param axis A unit vector defining the axis about which to rotate by angle theta
         * @param theta The angle at which to rotate about vector given by axis
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
     * \brief Class to hold information about the Backbone
     */
    class Backbone {
    public:

        Backbone() : backbone{}, file_path{}, interconnects{}, linker{}, vector_atom_deleted{} {}

        /**
         * \brief Backbone unit
         * @param file_path The path to the file containing the backbone molecule
         * @param interconnects The atom indices that define the periodic conditions between backbones
         * @param linker The atom indices used to align and connect backbone to base in the form
         * {atom to which base is bonded, hydrogen used to define vector for alignment to be deleted}
         */
        Backbone(std::string file_path, std::array<unsigned, 2> interconnects, std::array<unsigned,2> linker); 
            

        /**
         * \brief Gives the pointer to an atom that is the head from interconnects{head, tail}
         * @return The atom pointer from the backbone OBMol object
         */
        OpenBabel::OBAtom* getHead() {
            return backbone.GetAtom(interconnects[0]);
        }

        /**
         * \brief Gives the pointer to an atom that is the tail from interconnects{head, tail}
         * @return Pointer to the atom for the tail
         */
        OpenBabel::OBAtom* getTail() {
            return backbone.GetAtom(interconnects[1]);
        }

        /**
         * \brief Get the linker atom pointer
         * @return The atom that is the one linking to the Base
         */
        OpenBabel::OBAtom* getLinker() {
            return backbone.GetAtom(linker[0]);
        }

        /**
         * \brief Get the vector atom pointer (which is probably a hydrogen)
         * @return The atom that defines the vector from the
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
         * \brief Centers the molecule. Basically just an alias of the Center() function from OpenBabel.
         */
        void center() {
            backbone.Center();
        }

        /**
         * \brief Rotates the molecule by a matrix. Basically just an alias of the Rotate() function from OpenBabel.
         * @param rot The matrix by which to rotate the molecule
         */
        void rotate(double* rot) {
            backbone.Rotate(rot);
        }

        /**
         * \brief Translates the molecule by a vector. Basically just an alias of the Translate() function from OpenBabel.
         * @param vec The vector by which to translate the molecule
         */
        void translate(OpenBabel::vector3 vec) {
            backbone.Translate(vec);
        }

        /**
         * \brief Gives a copy of the molecule in the Backbone
         * @return A copy of the backbone molecule
         */
        OpenBabel::OBMol getMolecule() {
            return backbone;
        }

        /**
         * \brief Deletes the atom from \code{getVector()} safely. If the atom is already deleted, nothing happens.
         */
        void deleteVectorAtom();

       /**
        * \brief Does some basic sanity checks (such as whether or not the indices of the atom are within the range
        * of the molecule
        */
        void validate();

        std::array<unsigned , 2> interconnects,   //!< \brief { head, tail }
                linker,                             //!< \brief { atom index connecting to backbone, hydrogen defining vector }
                new_interconnects;                  //!< \brief Fixed Bonds
        OpenBabel::OBMol backbone;                  //!< \brief The molecule for the backbone
        std::string file_path;
        bool vector_atom_deleted;                   //!< \brief Whether or not the atom from \code{getVector()} has been deleted
    };

    class Base;

    /**
     * \brief Class to fully define bases (i.e. Adenine, Cytosine)
     */
    class Base {

    public:

        Base() : name{}, code{}, linker{}, base{}, vector_atom_deleted{}, pair_name{} {};

        /**
         * \brief Create Base from basic set of parameters
         * @param nameT The name of the base, full name
         * @param codeT The code-name of the base
         * @param file_path The path to the file containing the base molecule
         * @param linkerT An array that contains how the base is to link to a backbone
         * @param pairT The name of the pair from other bases which is empty by default
         * in the form {base linker, hydrogen atom defining vector}
         */

        Base(std::string nameT, std::string codeT, std::string file_path, std::array<std::size_t, 2> linkerT,
             std::string pairT = "");

        /**
         * \brief Gives the atom of the base that connects directly to the backbone
         * @return A pointer to the atom that connects
         */
        OpenBabel::OBAtom* getLinker() {
            return base.GetAtom(static_cast<unsigned>(linker[0]));
        }

        /**
         * \brief Gives the (most likely hydrogen) atom of the base connected to the atom from \code{getLinker()} which defines how
         * the base connects
         * @return
         */
        OpenBabel::OBAtom* getVector() {
            if (!vector_atom_deleted)
                return base.GetAtom(static_cast<unsigned>(linker[1]));
            else
                return nullptr;
        }

        /**
         * \brief Returns a copy of the base molecule
         * @return A copy of the base molecule
         */
        OpenBabel::OBMol getMolecule() {
            return base;
        }

        /**
         * \brief Gives the code of the Base
         * @return The code of the base
         */
        std::string getCode() {
            return code;
        }

        /**
         * \brief Gives the full name of the base
         * @return The full name of the base
         */
        std::string getName() {
            return name;
        }

        /**
         * \brief Deletes the atom from \code{getVector()} safely. If the atom is already deleted, nothing happens.
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
         * \brief Get the name of the pair base
         */
        std::string getBasePairName() {
            return pair_name;
        }

        /**
         * \brief Does some basic sanity checks (such as whether or not the indices of the atom are within the range
         * of the molecule
         */
        void validate();

        std::string name,                               //!< \brief Full name of base (i.e. Adenine)
                    code,                               //!< \brief Three character code to define base (Adenine: ADE)
                    pair_name,                          //!< \brief Name of the pair base
                    file_path;                          //!< \brief Path to a file containing the base
        OpenBabel::OBMol base;                          //!< \brief Pointer to OBMol defining the base
        std::array<std::size_t, 2 > linker;             //!< \brief Holds indices for atoms forming a vector to connect to backbone {linker, hydrogen}
        bool vector_atom_deleted;                       //!< \brief Whether or not the \code{getVector()} atom was deleted
    };

    /**
     * \brief A basic class that only contains a vector of Bases
     */
    struct Bases {
    public:
        /**
         * \brief Basic constructor for the Bases so we can feed it the base parameters created in python
         * @param py_bases The vector containing the information about the Bases needed
         */
        Bases(std::vector<Base> py_bases);

        Bases() {};

        PNAB::Base getBaseFromName(std::string name) {
            for (auto v : bases) {
                if (v.getName().find(name) != std::string::npos)
                    return v;
            }
            std::cerr << "Base \"" << name << "\" does not exists in list of bases. Please check input file."
                      << std::endl;
            throw 1;
        }

        std::vector<Base> getBasesFromStrand(std::vector<std::string> strand);
        std::vector<Base> getComplimentBasesFromStrand(std::vector<std::string> strand);

    private:
        std::vector<Base> bases;                        //!< \brief The vector of bases
        bool all_bases_pair;
        std::map<std::string, PNAB::Base> name_base_map;
    };

    /**
     * \brief Class to hold bases with backbones attached, along with associated necessary information
     */
    class BaseUnit {
    public:
        BaseUnit(Base b, Backbone backbone);
        const OpenBabel::OBMol getMol() {
            return unit;
        }

        const std::array< std::size_t, 2 > getBaseIndexRange() {
            return base_index_range;
        };

        const std::array< std::size_t, 2 > getBaseBondIndices() {
            return base_bond_indices;
        }

        const std::array< std::size_t, 2 > getBackboneIndexRange() {
            return backbone_index_range;
        };

        const std::array< std::size_t, 2 > getBackboneLinkers() {
            return backbone_interconnects;
        };

        std::size_t getBaseConnectIndex() {
            return base_connect_index;
        }



    private:
        OpenBabel::OBMol unit;                                                       //!< \brief Holds molecule containing base with backbone attached
        OpenBabel::OBMol base_;
        std::array< std::size_t, 2 > base_index_range,                               //!< \brief Range of indices of the unit that are a part of the base, [start, stop]
                                     backbone_index_range;                           //!< \brief Range of indices of the unit that are a part of the backbone, [start, stop]
        std::size_t base_connect_index;                                              //!< \brief Atom index where Backbone connects to Base (the Base atom)
        std::array< std::size_t, 2 > backbone_interconnects;                         //!< \brief Atom indices defining where backbone connects
        std::array< std::size_t, 2 > base_bond_indices;                              //!< \brief Index range corresponding to bonds in the base of the BaseUnit
    };

    /**
     * \brief Class to contain important information for an individual conformer, specifically includes detailed
     * information about the energy components important for distinguishing between different conformers.
     * All energy values are a result of the total of each source of energy divided by the number of \code{UnitChain}s
     * used during testing to give a comparable value.
     */
    struct ConformerData {
        double *coords,                             //!< \brief Pointer to array containing coordinates of all atoms in molecule chain
               *monomer_coord,                      //!< \brief Pointer to array containing coordinates of a single monomer
                distance,                           //!< \brief distance between \code{interconnects} in \code{Backbone} for successive \code{UnitChain}s
                bondE,                              //!< \brief Total energy of bonds
                angleE,                             //!< \brief Total energy of angles of the conformation divided by number of \code{UnitChain}s
                VDWE,                               //!< \brief Total van Der Wals Energy
                total_energy,                       //!< \brief Total energy of the conformation divided by number of \code{UnitChain}s in chain tested
                rmsd;                               //!< \brief Root-mean square distance from lowest energy conformer
        std::size_t index;                          //!< \brief The index of the conformer
        bool chain_coords_present;                  //!< \brief Have the chain coordinates in coord been deleted?

        /**
         * \brief Used for simple sorting based on total energy of the conformer
         * @param cd The ConformerData element to compare current element to
         * @return True if the other ConformerData has greater total energy, false otherwise
         */
        bool operator<(const ConformerData &cd) const {
            return (total_energy < cd.total_energy);
        }
    };

    const auto DEGSS = 0.017453292519943295769236907684886127134428718885417254560971914401710;
}

#endif //PNAB_CONTAINERS_H
