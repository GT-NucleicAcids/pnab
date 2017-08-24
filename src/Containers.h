//
// Created by jbarnett8 on 8/22/17.
//

#ifndef PNAB_CONTAINERS_H
#define PNAB_CONTAINERS_H

#include <string>
#include <openbabel/mol.h>

class RuntimeParameters {

private:
    // Helical parameters
    double  H_Rise,
            mZ[9],
            mY[9],
            mX[9];
    std::array< double, 3 > angles;             // { inclination (x), tip (y), h_twist (z) }
    std::array< double, 2 > displacement;       // { x displacement, y displacement }

    // Backbone parameters
    std::array< std::size_t, 2 > interconnects, // { head, tail }
            baseconnect,                        // { connect, vector }
            new_interconnects;
    OpenBabel::OBMol *mol;	                    // Fixed Bonds

    // Energy parameters
    std::vector< double > energy_filter;        // { max total E, max angle E, max bond E, max VDW E, max Torsion E }
    double max_distance;

    // Force Field Parameters
    std::string type,
                parameter_file;
    double base_to_backbone_bond_length;

    // Search algorithm
    std::size_t num_steps,
                dihedral_discretization,		// only for weighted methods
                angleStepSize,
                chain_length,
    std::string algorithm;
};

struct base { 		                            // Class to fully define bases (i.e. Adenine, Cytosine)
    std::string  name,		                    // Full name of base (i.e. Adenine)
                 code,		                    // Three character code to define base (Adenine: ADE)
                 pair;		                    // String code for the pair (i.e. would be Thymine for Adenine)
    OpenBabel::OBMol *mol;			            // Pointer to OBMol defining the base
    std::size_t symmetry, 	                    // Discrete rotational symmetry order (mirror symmetry: n = 2)
                connect,				        // Index of the atom that connects to backbone
                vector, 				        // Index of hydrogen connected to connect, defines vector
                numAtoms;				        // Number of atoms in base to help offset indices of backbone atoms
};

struct conformerData {
    double *coords,
            distance,
            energy,
            angleE,
            bondE,
            VDWE,
            torsionE,
            totTorsionE,
            rmsd;
    int index;
    bool accept;
    bool operator < (const conformerData& cd) const {
        return (energy < cd.energy);
    }
};

#endif //PNAB_CONTAINERS_H
