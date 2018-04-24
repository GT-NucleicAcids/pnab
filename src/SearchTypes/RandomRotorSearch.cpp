//
// Created by jbarnett8 on 9/17/17.
//

#include "RandomRotorSearch.h"

using namespace OpenBabel;
using namespace PNAB;
using namespace std;

RandomRotorSearch::RandomRotorSearch(PNAB::RuntimeParameters &runtime_params_v, UnitChain &uc_v) {
//    runtime_params = runtime_params_v;
//    uc = uc_v;
//    chain = uc.getChainMol();
//    backbone_interlinks = uc.getBackboneLinkerIndices();
//    max_backbone_interlink_distance = runtime_params.max_distance;
//    rng.seed(std::random_device()());
//    range_vec = {runtime_params.x_disp, runtime_params.y_disp, runtime_params.rise,
//                 runtime_params.inclination, runtime_params.tip, runtime_params.twist};
//    translation = {range_vec[0].d[0], range_vec[1].d[0], range_vec[2].d[0]};
//    rotation = {range_vec[3].d[0], range_vec[4].d[0], range_vec[5].d[0]};
//    //vector<double&> val_vec = {translation[0], translation[1], translation[2], rotation[0], rotation[1], rotation[2]};
//    for (auto i : rotation)
//        cout << i << ", ";
//    cout << endl;
//    size_t num_base_units = backbone_interlinks.size();
//    double cs = cos(num_base_units * rotation[2] * DEG_TO_RAD),
//           sn = sin(num_base_units * rotation[2] * DEG_TO_RAD),
//           m_uc_twist[3][3] = {{cs, -sn, 0}, {sn, cs, 0}, {0, 0, 1}};
//    uc_twist = matrix3x3(m_uc_twist);
//    uc_rise = vector3(0, 0, num_base_units * translation[2]);

//    for (auto v : range_vec) {
//        if (v.d.size() >= 2 && v.type == RuntimeParameters::RANGE) {
//            uniform_real_distribution<double> dist(v.d[0], v.d[1]);
//            num_range_pair_vec.push_back(
//                    make_pair<double&, std::uniform_real_distribution<double>>())
//            dist_vec.push_back(dist);
//        } else {
//            uniform_real_distribution<double> dist;
//            dist_vec.push_back(dist);
//        }
//    }


}

OBMol RandomRotorSearch::search() {
//    double *coords = uc.getCoordinates();
//    uniform_real_distribution<double> dist(0, 2 * M_PI);
//
//    OBBitVec fix_atoms(chain.NumAtoms());
//    OBBitVec fix_bonds(chain.NumBonds());
//    auto base_indices = uc.getBaseIndices();
//    auto bond_indices = uc.getBaseBondIndices();
//    OBRotorList rl;
//    for (auto v : base_indices)
//        for (unsigned i = static_cast<unsigned>(v[0]); i <= v[1]; ++i)
//            fix_atoms.SetBitOn(i);
//    for (auto v : bond_indices)
//        for (unsigned i = static_cast<unsigned >(v[0]); i <= v[1]; ++i)
//            fix_bonds.SetBitOn(i);
//    rl.Setup(chain);
//    rl.SetFixAtoms(fix_atoms);
//    //rl.SetFixedBonds(fix_bonds);
//    rl.SetRotAtomsByFix(chain);
//    //rl.Setup(chain);
//
//    size_t search_size = runtime_params.num_steps;
//    OBRotorIterator ri;
//    for (size_t search_index = 0; search_index < search_size; ++search_index) {
//        OBRotor *r = rl.BeginRotor(ri);
//        for (int i = 0; i < rl.Size(); i++, r = rl.NextRotor(ri)) {
//            r->SetToAngle(coords, dist(rng));
//        }
//        if (passesBackboneLinkerDistanceFilter(coords)) {
//            vector<double> coord_vec(3*chain.NumAtoms());
//            memcpy(coord_vec.data(), coords, 3*chain.NumAtoms()*sizeof(double));
//            coords_vec.push_back(coord_vec.data());
//        }
//        if (search_index % 100000 == 0)
//            cout << 100 * static_cast<double>(search_index) / search_size << "%\t\tAccepted: " << coords_vec.size() << endl;
//    }
//    if (coords_vec.size() > 0) {
//        chain.SetCoordinates(coords_vec[0]);
//    } else {
//        chain.SetCoordinates(coords);
//    }
//    return chain;
}


bool RandomRotorSearch::passesBackboneLinkerDistanceFilter(double *coord) {
    chain.SetCoordinates(coord);
    if (backbone_interlinks.size() >= 2)
        for (unsigned i = 1; i < backbone_interlinks.size(); ++i) {
            OBAtom *head = chain.GetAtom(backbone_interlinks[i - 1][1]);
            if (head->GetDistance(backbone_interlinks[i][0]) > max_backbone_interlink_distance)
                return false;
        }
    vector3 tail = chain.GetAtom(backbone_interlinks[0][0])->GetVector();
    vector3 head = chain.GetAtom(backbone_interlinks.back()[1])->GetVector();

    // cout << "tail: " << tail << ", head: " << head << endl;
    // cout << "\tDistance: " << head.distSq(tail) << ", square-root: " << sqrt(head.distSq(tail)) << endl;
    tail *= uc_twist;
    tail += uc_rise;
    // cout << "Distance: " << head.distSq(tail) << ", square-root: " << sqrt(head.distSq(tail)) << endl;
    if (sqrt(head.distSq(tail)) > max_backbone_interlink_distance)
        return false;
    return true;
}

void RandomRotorSearch::updateGeometry() {
//    for (unsigned i = 0; i < range_vec.size(); ++i) {
//        auto v = range_vec[i];
//        if (v.type == RuntimeParameters::RANGE)
//            if (i > 2)
//                rotation[i] = dist_vec[i](rng);
//            else
//                translation[i] = dist_vec[i](rng);
//    }
//    uc.updateHelicalParameters(translation, rotation);
}

//void RandomRotorSearch::addRange(RuntimeParameters::RangeField &field, double &val) {
//    if (field.d.size() >= 2 && field.type == RuntimeParameters::RANGE) {
//        uniform_real_distribution<double> dist(field.d[0], field.d[1]);
//        //auto& rng_tmp = rng;
//        std::mt19937_64 rng_tmp;
//        rng_tmp.seed(std::random_device()());
//        //auto a = dist(rng_tmp);
//        params_functions.push_back([val, dist, &rng_tmp](){ val = dist(rng_tmp); });
//    }
//}
