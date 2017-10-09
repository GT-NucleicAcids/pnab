//
// Created by jbarnett8 on 9/4/17.
//

#include "UnitChain.h"
#include <openbabel/math/align.h>

using namespace std;
using namespace PNAB;
using namespace OpenBabel;

UnitChain::UnitChain(std::vector< std::string > chain_string, Bases bases, Backbone backbone) {
    for (auto chain_string_n : chain_string ) {
        stringstream ss(chain_string_n);
        std::vector< BaseUnit > base_units_n;
        string tok;
        while (getline(ss, tok, ',')) {
            bool base_found = false;
            transform(tok.begin(), tok.end(), tok.begin(), ::tolower);
            for (auto b : bases.bases) {
                if (tok.find(b.getCode()) != string::npos) {
                    BaseUnit bu(b, backbone);
                    base_units_n.push_back(bu);
                    base_found = true;
                    break;
                }
            }
            if (!base_found) {
                cerr << "Missing base with code \"" << tok << "\". Please check base codes "
                     << "and string of base codes provided in input file. " << endl;
                exit(1);
            }
        }
        base_units.push_back( base_units_n );
        if (base_units.size() > 0 && base_units[0].size() != base_units_n.size()) {
            cerr << "Length of chain string \"" << chain_string_n << "\" does not match "
                 << "length of first chain string (" << base_units[0].size() << " base units). "
                 << "The dimensions must match." << endl;
            exit(1);
        }
    }
}