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
    unsigned symmetry, 	            // Discrete rotational symmetry order (mirror symmetry: n = 2)
             connect,				// Index of the atom that connects to backbone
             vector, 				// Index of hydrogen connected to connect, defines vector
             numAtoms;
};

#endif //PNAB_CONTAINERS_H
