//
// Created by jbarnett8 on 8/22/17.
//

#ifndef PNAB_CONTAINERS_H
#define PNAB_CONTAINERS_H

#include <string>
#include <openbabel/mol.h>
#include "FileParser.h"

namespace PNAB {

    /**
     * \brief A class for holding all necessary parameters for conformational searches
     */
    class RuntimeParameters {

    public:
        /**
         * \brief Empty constructor. Should not generally be used.
         */
        RuntimeParameters() : mZ{}, mY{}, mX{}, inclination{}, tip{}, twist{}, x_disp{}, y_disp{},
                              energy_filter{}, max_distance(0), type(), parameter_file(),
                              base_to_backbone_bond_length(-1), num_steps(0),
                              dihedral_discretization(0), angleStepSize(0), chain_length(3),
                              algorithm() {};

        /**
         * \brief Constructs RuntimeParameters from a FileParser object that is already constructed from an input file
         * @param sp The FileParser that contains all of the desired information
         */
        RuntimeParameters(FileParser &sp);

        /**
         * \brief Basic sanity checks for values. Exits the program if errors are encountered.
         */
        void validate();

    private:
        // Helical parameters
        double mZ[9],                       //!< Holds rotation matrix about x-axis
                mY[9],                      //!< Holds rotation matrix about y-axis
                mX[9];                      //!< Holds rotation matrix about z-axis

        /**
         * \brief Holds values that tell us whether the geomtric parameters in range_field are defined as a single value or a range
         * which would imply we must also search over this range
         */
        enum RangeType {
            SINGLE,                             //!< If only a single value is defined
            RANGE                               //!< If a range of values is defined (begin, end)
        };

        /**
         * \brief A struct for holding geometric parameters as well as the state (Single or Range).
         */
        struct RangeField {
            std::vector<double> d;
            RangeType type;
        };

        RangeField inclination,                 //!< Rotation about x-axis (applied once to BaseUnit)
                    tip,                        //!< Rotation about y-axis (applied once to BaseUnit)
                    twist,                      //!< Rotation about z-axis (applied n times to BaseUnit, where n is number of bases in BaseUnit)
                    rise,                       //!< Displacement in z direction (applied n times to BaseUnit)
                    x_disp,                     //!< Displacement in x direction (applied once to BaseUnit)
                    y_disp;                     //!< Displacement in y direction (applied once to BaseUnit)

        // Energy parameters
        std::vector<double> energy_filter;      //!< { max total E, max angle E, max bond E, max VDW E, max Torsion E }
        double max_distance;                    //!< The maximum distance between head and tail of successive UnitChains that is accepted

        // Force Field Parameters
        std::string type,                       //!< The type of the ForceField such as "GAFF" or "MMFF94"
                parameter_file;                 //!< An additional parameter file in case missing bonds, angles, torsions, etc. are missing
        double base_to_backbone_bond_length;    //!< The bond length between Base and Backbone, can be left to default which is original distance of the distance between the atom of getLinker() and getVector() of the Backbone

        // Search algorithm
        std::size_t num_steps,                  //!< Determines how many points are sampled in the Monte Carlo searches
                dihedral_discretization,        //!< The number of points to break up the total number of dihedral angle (step becomes 360 deg / dihedral_discretization)
                angleStepSize,                  //!< The step size for a systematic search
                chain_length;                   //!< The number of UnitChains to use for energy testing
        std::string algorithm;                  //!< The algorithm to use.

        /**
         * \brief Checks to see of the vector is a range for a single value then sets the parameters inside of the
         * the \code{strct}.
         * @param strct Contains the data for the range and a range type
         */
        void checkRangeField(RangeField &strct) {
            const std::size_t max_size = 2;
            std::size_t d_size = strct.d.size();
            if (d_size > max_size) {
                std::cerr << "Vector size is too large" << std::endl;
                exit(1);
            }
            if (d_size == 1)
                strct.type = SINGLE;
            else
                strct.type = RANGE;
        }
    };

    /**
     * \brief Class to hold information about the Backbone
     */
    class Backbone {
    public:

        /**
         * \brief Backbone unit
         * @param molP The molecule containing the molecule
         * @param interconnects The atom indices that define the periodic conditions between backbones
         * @param linkerT The atom indices used to align and connect backbone to base in the form
         * {atom to which base is bonded, hydrogen used to define vector for alignment to be deleted}
         */
        Backbone(OpenBabel::OBMol &molP, std::array<unsigned, 2> interconnects, std::array<unsigned, 2> linker) {
            backbone = OpenBabel::OBMol(molP);
            this->interconnects = interconnects;
            this->linker = linker;
            vector_atom_deleted = false;
        }

        /**
         * \brief Create a Backbone from from a FileParser
         * @param fp The FileParser containing all the information for the Backbone
         */
        Backbone(FileParser &fp);


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
            else
                return nullptr;
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

    private:
        std::array<unsigned , 2> interconnects,   //!< { head, tail }
                linker,                             //!< { atom index connecting to backbone, hydrogen defining vector }
                new_interconnects;                  //!< Fixed Bonds
        OpenBabel::OBMol backbone;                  //!< The molecule for the backbone
        bool vector_atom_deleted;                   //!< Whether or not the atom from \code{getVector()} has been deleted
    };

    /**
     * \brief Class to fully define bases (i.e. Adenine, Cytosine)
     */
    class Base {

    public:
        /**
         * \brief Create Base from basic set of parameters
         * @param nameT The name of the base, full name
         * @param codeT The code-name of the base
         * @param molT The molecule containing the base's chemical makeup
         * @param linkerT An array that contains how the base is to link to a backbone
         * in the form {base linker, hydrogen atom defining vector}
         */
        Base(std::string nameT, std::string codeT, OpenBabel::OBMol &molT, std::array<std::size_t, 2> linkerT) {
            name = nameT;
            code = codeT;
            linker = linkerT;
            base = OpenBabel::OBMol(molT);
            vector_atom_deleted = false;
            validate();
        }

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
         * \brief Does some basic sanity checks (such as whether or not the indices of the atom are within the range
         * of the molecule
         */
        void validate();

    private:
        std::string name,                               //!< Full name of base (i.e. Adenine)
                code;                                   //!< Three character code to define base (Adenine: ADE)
        OpenBabel::OBMol base;                          //!< Pointer to OBMol defining the base
        std::array<std::size_t, 2 > linker;             //!< Holds indices for atoms forming a vector to connect to backbone {linker, hydrogen}
        bool vector_atom_deleted;                       //!< Whether or not the \code{getVector()} atom was deleted
    };

    /**
     * \brief A basic class that only contains a vector of Bases
     */
    struct Bases {
    public:
        /**
         * \brief Basic constructor for the Bases so we can feed it the FileParser as an input to make the Base
         * @param fp The FileParser that contains all of the information about the Bases needed
         */
        Bases(FileParser &fp);

        std::vector<Base> bases;                        //!< The vector of bases
    };

    /**
     * \brief Class to hold bases with backbones attached, along with associated necessary information
     */
    class BaseUnit {
    public:
        BaseUnit(Base b, Backbone backbone);
        const OpenBabel::OBMol* getMol() {
            return &unit;
        }

    private:
        OpenBabel::OBMol unit;                      //!< Holds molecule containing base with backbone attached
        OpenBabel::OBAtom *baseAtomBegin;
        OpenBabel::OBAtom *baseAtomEnd;
        OpenBabel::OBAtom *backboneAtomBegin;
        OpenBabel::OBAtom *backboneAtomEnd;
    };

    /**
     * \brief Class to contain important information for an individual conformer, specifically includes detailed
     * information about the energy components important for distinguishing between different conformers.
     * All energy values are a result of the total of each source of energy divided by the number of \code{UnitChain}s
     * used during testing to give a comparable value.
     */
    struct ConformerData {
        double *coords,                             //!< Pointer to array containing coordinates of all atoms in molecule
                distance,                           //!< distance between \code{interconnects} in \code{Backbone} for successive \code{UnitChain}s
                total_energy,                       //!< Total energy of the conformation divided by number of \code{UnitChain}s in chain tested
                angleE,                             //!< Total energy of angles of the conformation divided by number of \code{UnitChain}s
                bondE,                              //!< Total energy of bonds
                VDWE,                               //!< Total van Der Wals Energy
                torsionE,                           //!< Torsion energy due to backbone
                totTorsionE,                        //!< Total torsional energy
                rmsd;                               //!< Root-mean square distance from lowest energy conformer
        int index;                                  //!< The index of the conformer
        bool accept;                                //!< Did the conformer pass the energy filter

        /**
         * \brief Used for simple sorting based on total energy of the conformer
         * @param cd The ConformerData element to compare current element to
         * @return True if the other ConformerData has greater total energy, false otherwise
         */
        bool operator<(const ConformerData &cd) const {
            return (total_energy < cd.total_energy);
        }
    };

}

#endif //PNAB_CONTAINERS_H
