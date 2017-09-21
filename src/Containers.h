//
// Created by jbarnett8 on 8/22/17.
//

#ifndef PNAB_CONTAINERS_H
#define PNAB_CONTAINERS_H

#include <string>
#include <openbabel/mol.h>

namespace PNAB {

    class RuntimeParameters {

    private:
        // Helical parameters
        double H_Rise,
                mZ[9],
                mY[9],
                mX[9];
        std::array<double, 3> angles;             // { inclination (x), tip (y), h_twist (z) }
        std::array<double, 2> displacement;       // { x displacement, y displacement }

        // Energy parameters
        std::vector<double> energy_filter;        // { max total E, max angle E, max bond E, max VDW E, max Torsion E }
        double max_distance;

        // Force Field Parameters
        std::string type,
                parameter_file,
                algorithm;
        double base_to_backbone_bond_length;

        // Search algorithm
        std::size_t num_steps,
                dihedral_discretization,        // only for weighted methods
                angleStepSize,
                chain_length;
    };

    /**
     * Class to hold information about the Backbone
     */
    class Backbone {
    public:

        Backbone(OpenBabel::OBMol &molP, std::array<unsigned, 2> interconnectsT, std::array<unsigned, 2> linkerT);

        OpenBabel::OBAtom* getHead();

        OpenBabel::OBAtom* getTail();

        OpenBabel::OBAtom* getLinker();

        OpenBabel::OBAtom* getVector();

        void center();
        void rotate(double* rot);
        void translate(OpenBabel::vector3 vec);

        OpenBabel::OBMol getMolecule();

        void deleteVectorAtom();

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
        Base(std::string nameT, std::string codeT, OpenBabel::OBMol &molT, std::array<unsigned, 2> linkerT);

        OpenBabel::OBAtom* getLinker();

        OpenBabel::OBAtom* getVector();

        OpenBabel::OBMol getMolecule();

        std::string getCode();

        std::string getName();

        void deleteVectorAtom();

    private:
        std::string name,                               //!< Full name of base (i.e. Adenine)
                code;                                   //!< Three character code to define base (Adenine: ADE)
        OpenBabel::OBMol base;                          //!< Pointer to OBMol defining the base
        std::array< unsigned, 2 > linker;               //!< Holds indices for atoms forming a vector to connect to backbone {linker, hydrogen}
        bool vector_atom_deleted;
    };

    /**
     * Class to hold bases with backbones attached, along with associated necessary information
     */
    class BaseUnit {
    public:
        BaseUnit(Base b, Backbone bb);

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
