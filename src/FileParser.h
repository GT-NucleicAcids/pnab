//
// Created by jbarnett8 on 9/20/17.
//

#ifndef PNAB_FILEPARSER_H
#define PNAB_FILEPARSER_H

#include <string>
#include <vector>

class FileParser {

public:
    FileParser();

    void registerCategory();
    void registerField();

private:

};

class Category {
    std::string name;
    std::vector< Field > fields;
};

class Field {
    std::string name;
};

#endif //PNAB_FILEPARSER_H
