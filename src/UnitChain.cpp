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
        transform(chain_string_n.begin(), chain_string_n.end(), chain_string_n.begin(), ::toupper);
        stringstream ss(chain_string_n);
        std::vector< BaseUnit > base_units_n;
        string tok;
        while (getline(ss, tok, ',')) {
            bool base_found = false;
            for (auto b : bases.bases) {
                if (tok.find(b.getCode()) != string::npos) {
                    BaseUnit bu(b, backbone);
                    base_units_n.push_back(bu);
                    addToChain(bu);
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

void UnitChain::addToChain(BaseUnit &bu) {
    // We need the number of atoms already present to offset absolute value of indices
    size_t num_atoms_before = chain.NumAtoms();
    size_t num_bonds_before = chain.NumBonds();

    // Get the Base atom ranges then add the offset, then store
    auto arr = bu.getBaseIndexRange();
    for (auto &v : arr)
        v += num_atoms_before;
    base_indices.push_back(arr);

    // Get the Backbone ranges then add the offset, then store
    arr = bu.getBackboneIndexRange();
    for (auto &v : arr)
        v += num_atoms_before;
    backbone_indices.push_back(arr);

    // The the inter-Backbone linker indices and add the offset, then store
    arr = bu.getBackboneLinkers();
    for (auto &v : arr)
        v+= num_atoms_before;
    backbone_linker_indices.push_back(arr);

    // The base bond indices and the offset, then store
    arr = bu.getBaseBondIndices();
    for (auto &v : arr)
        v+= num_bonds_before;
    base_bond_indices.push_back(arr);

    // Add the molecule
    chain += bu.getMol();

    size_t num_atoms = chain.NumAtoms();

    // Base is added first, so it is the first index of the BaseUnit while the last atom in Backbone_indices unit is the
    // last atom in the BaseUnit
    base_unit_indices.push_back({base_indices.back()[0], backbone_indices.back()[1]});

    // Store the updated original coordinates
    double *coords = chain.GetCoordinates();
    original_coords.resize(3 * num_atoms);
    memcpy(original_coords.data(), coords, 3 * num_atoms * sizeof(double));
}

void UnitChain::updateHelicalParameters(std::array<double, 3> &translation, std::array<double, 3> &rotation) {
    chain.SetCoordinates(original_coords.data());
    array< double, 3 > deg_rot = {DEGSS * rotation[0], DEGSS * rotation[1], DEGSS * rotation[2]};
    array< double, 9 > rot_mat = {1,               0,                0,
                                  0, cos(deg_rot[0]), -sin(deg_rot[0]),
                                  0, sin(deg_rot[0]),  cos(deg_rot[0])};
    chain.Rotate(rot_mat.data());
    rot_mat = {cos(deg_rot[1]), 0, -sin(deg_rot[1]),
                             0, 1,                0,
               sin(deg_rot[1]), 0,  cos(deg_rot[1])};
    chain.Rotate(rot_mat.data());
    chain.Translate(translation.data());
    double *coords = chain.GetCoordinates();
    int unit_count = -1;
    for (auto range : base_unit_indices) {
        // cout << "range[0], range[1]: " << range[0] << ", " << range[1] << endl;
        for (size_t i = range[0] - 1; i <= range[1] - 1; ++i) {
            rot_mat = {cos(unit_count * deg_rot[2]), -sin(unit_count * deg_rot[2]),  0,
                       sin(unit_count * deg_rot[2]),  cos(unit_count * deg_rot[2]),  0,
                                                  0,                             0,  1};
            double xT = coords[3 * i],
                   yT = coords[3 * i + 1],
                   zT = coords[3 * i + 2];
            double x  = rot_mat[0] * xT + rot_mat[1] * yT + rot_mat[2] * zT,
                   y  = rot_mat[3] * xT + rot_mat[4] * yT + rot_mat[5] * zT,
                   z  = rot_mat[6] * xT + rot_mat[7] * yT + rot_mat[8] * zT;
            z += unit_count * translation[2];
            coords[3 * i    ] = x;
            coords[3 * i + 1] = y;
            coords[3 * i + 2] = z;
            // cout << coords[3 * i + 2] << ", ";
        }
        // cout << endl;
        ++unit_count;
    }
}

