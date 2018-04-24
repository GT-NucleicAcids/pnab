//
// Created by jbarnett8 on 8/24/17.
//

#ifndef PNAB_SEARCHTEMPLATE_H
#define PNAB_SEARCHTEMPLATE_H

#include "../FileParser.h"
#include "../Containers.h"
#include "../UnitChain.h"

/**
 * \brief Template class for all search algorithms
 *
 * This is the template class for all search algorithms defined for a unified interface. It makes it easier to use.
 */
class SearchAlgorithm {

public:
    void Setup(FileParser fp, PNAB::UnitChain uc);

};

#endif //PNAB_SEARCHTEMPLATE_H
