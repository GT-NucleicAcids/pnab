//
// Created by jbarnett8 on 8/22/17.
//

#ifndef PNAB_CONTAINERS_H
#define PNAB_CONTAINERS_H

#include <string>
#include <openbabel/mol.h>

class Base {

public:
    Base();
    ~Base();

private:
    std::string name,		        // Full name of base (i.e. Adenine)
                code;		        // Three character code to define base (Adenine: ADE)
    OpenBabel::OBMol *mol;			// Pointer to OBMol defining the base
    std::size_t symmetry, 	        // Discrete rotational symmetry order (mirror symmetry: n = 2)
                connect,			// Index of the atom that connects to backbone
                vector, 			// Index of hydrogen connected to connect, defines vector
                num_atoms;
};

class BasePair {

public:
    BasePair();
    ~BasePair();

private:
    std::vector< Base > bases;
};

class HelicalData {

public:
    HelicalData();
    ~HelicalData();

private:
    double  h_rise,
            inclination,
            tip,
            h_twist,
            x_disp,
            y_disp,
            mZ[9],
            mX[9],
            mY[9];


};

class BackboneParams {
    unsigned head,
             tail,
             connect,
             vector,
             newHead,
             newTail;
    OpenBabel::OBMol *mol;

};

class ConformerData {

public:
    ConformerData();
    ~ConformerData();

private:

    double *coords,
            distance,
            tot_energy,
            angle_energy,
            bond_energy,
            vdw_energy,
            torsion_energy,
            tot_torsion_energy,
            rmsd;
    std::size_t index;
    bool accept;
    bool operator < (const ConformerData& cd) const {
        return (tot_energy < cd.tot_energy);
    }
};

class ConformerFilterData {

public:
    ConformerFilterData();
    ~ConformerFilterData();

private:
    double max_energy,
           max_angle_energy,
           max_bond_energy,
           max_vdw_energy,
           max_torsion_energy,
           max_distance;
};

class ForceFieldData {

public:
    ForceFieldData();
    ~ForceFieldData();

private:
    std::string type,
                parameter_file;
    double base_to_backbone_bond_length;
};

class ConformerSearchData {
    std::size_t num_steps,
                num_dihedral_steps,
                angle_step_size,
                chain_length;
    std::string algorithm;
};

#endif //PNAB_CONTAINERS_H
