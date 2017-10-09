//
// Created by jbarnett8 on 8/22/17.
//

#ifndef PNAB_CONTAINERS_H
#define PNAB_CONTAINERS_H

#include <string>
#include <openbabel/mol.h>
#include "FileParser.h"

namespace PNAB {

    class RuntimeParameters {

    public:
        RuntimeParameters() : mZ{}, mY{}, mX{}, inclination{}, tip{}, twist{}, x_disp{}, y_disp{},
                              energy_filter{}, max_distance(0), type(), parameter_file(),
                              base_to_backbone_bond_length(-1), num_steps(0),
                              dihedral_discretization(0), angleStepSize(0), chain_length(3),
                              algorithm() {};

        RuntimeParameters(FileParser &sp);

        void validate();

    private:
        // Helical parameters
        double mZ[9],
                mY[9],
                mX[9];

        enum RangeType {
            INCL_SINGLE,
            INCL_RANGE,
            TIP_SINGLE,
            TIP_RANGE,
            TWIST_SINGLE,
            TWIST_RANGE,
            RISE_SINGLE,
            RISE_RANGE,
            X_DISP_SINGLE,
            X_DISP_RANGE,
            Y_DISP_SINGLE,
            Y_DISP_RANGE
        };

        struct range_field {
            std::vector<double> d;
            RangeType type;
        };

        range_field inclination,
                    tip,
                    twist,
                    rise,
                    x_disp,
                    y_disp;

        // Energy parameters
        std::vector<double> energy_filter;        // { max total E, max angle E, max bond E, max VDW E, max Torsion E }
        double max_distance;

        // Force Field Parameters
        std::string type,
                parameter_file;
        double base_to_backbone_bond_length;

        // Search algorithm
        std::size_t num_steps,
                dihedral_discretization,        // only for weighted methods
                angleStepSize,
                chain_length;
        std::string algorithm;

        /**
         * \brief Checks to see of the vector is a range for a single value then sets the parameters inside of the
         * the \code{strct}.
         * @param strct Contains the data for the range and a range type
         * @param single The type for when the strct contains only one value
         * @param range The type for when the strct contains two values
         */
        void checkRangeField(range_field &strct, RangeType single, RangeType range) {
            const std::size_t max_size = 2;
            std::size_t d_size = strct.d.size();
            if (d_size > max_size) {
                std::cerr << "Vector size is too large" << std::endl;
                exit(1);
            }
            if (d_size == 1)
                strct.type = single;
            else
                strct.type = range;
        }
    };

    /**
     * Class to hold information about the Backbone
     */
    class Backbone {
    public:

        /**
         * Backbone unit
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

        Backbone(FileParser &fp);


        /**
         * \brief Gives the pointer to an atom that is the head from interconnects{head, tail}
         * @return The atom pointer from the backbone OBMol object
         */
        OpenBabel::OBAtom* getHead() {
            return backbone.GetAtom(interconnects[0]);
        }

        OpenBabel::OBAtom* getTail() {
            return backbone.GetAtom(interconnects[1]);
        }

        OpenBabel::OBAtom* getLinker() {
            return backbone.GetAtom(linker[0]);
        }

        OpenBabel::OBAtom* getVector() {
            if (!vector_atom_deleted)
                return backbone.GetAtom(linker[1]);
            else
                return nullptr;
        }

        void center() {
            backbone.Center();
        }

        void rotate(double* rot) {
            backbone.Rotate(rot);
        }

        void translate(OpenBabel::vector3 vec) {
            backbone.Translate(vec);
        }

        OpenBabel::OBMol getMolecule() {
            return backbone;
        }

        void deleteVectorAtom();

        void validate();

    private:
        std::array<unsigned , 2> interconnects,   //!< { head, tail }
                linker,                             //!< { atom index connecting to backbone, hydrogen defining vector }
                new_interconnects;                  //!< Fixed Bonds
        OpenBabel::OBMol backbone;
        bool vector_atom_deleted;
    };

    /**
     * Class to fully define bases (i.e. Adenine, Cytosine)
     */
    class Base {

    public:
        Base(std::string nameT, std::string codeT, OpenBabel::OBMol &molT, std::array<std::size_t, 2> linkerT) {
            name = nameT;
            code = codeT;
            linker = linkerT;
            base = OpenBabel::OBMol(molT);
            vector_atom_deleted = false;
            validate();
        }

        OpenBabel::OBAtom* getLinker() {
            return base.GetAtom(static_cast<unsigned>(linker[0]));
        }

        OpenBabel::OBAtom* getVector() {
            if (!vector_atom_deleted)
                return base.GetAtom(static_cast<unsigned>(linker[1]));
            else
                return nullptr;
        }

        OpenBabel::OBMol getMolecule() {
            return base;
        }

        std::string getCode() {
            return code;
        }

        std::string getName() {
            return name;
        }

        void deleteVectorAtom() {
            if (!vector_atom_deleted) {
                size_t id = getLinker()->GetId();
                base.DeleteAtom(getVector());
                linker = {base.GetAtomById(id)->GetIdx(), 0};
                vector_atom_deleted = true;
            }
        }

        void validate();

    private:
        std::string name,                               //!< Full name of base (i.e. Adenine)
                code;                                   //!< Three character code to define base (Adenine: ADE)
        OpenBabel::OBMol base;                          //!< Pointer to OBMol defining the base
        std::array<std::size_t, 2 > linker;               //!< Holds indices for atoms forming a vector to connect to backbone {linker, hydrogen}
        bool vector_atom_deleted;
    };

    struct Bases {
    public:
        Bases(FileParser &fp);
        std::vector<Base> bases;
    };

    /**
     * Class to hold bases with backbones attached, along with associated necessary information
     */
    class BaseUnit {
    public:
        BaseUnit(Base b, Backbone bb);
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
     * Class to contain important information for an individual conformer, specifically includes detailed
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
                VDWE,
                torsionE,
                totTorsionE,
                rmsd;
        int index;
        bool accept;

        bool operator<(const ConformerData &cd) const {
            return (total_energy < cd.total_energy);
        }
    };

}

#endif //PNAB_CONTAINERS_H
