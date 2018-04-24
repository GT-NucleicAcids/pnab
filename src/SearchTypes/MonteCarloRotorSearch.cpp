//
// Created by jbarnett8 on 2/26/18.
//

#include "MonteCarloRotorSearch.h"

using namespace PNAB;
using namespace std;
using namespace OpenBabel;

MonteCarloRotorSearch::MonteCarloRotorSearch(RuntimeParameters &runtime_params, Base base_a,
                                             Base base_b, Backbone backbone,
                                             HelicalParameters &helical_params) {
    runtime_params_ = runtime_params;
    base_a_ = base_a;
    base_b_ = base_b;
    helical_params_ = helical_params;
    backbone_ = backbone;
    auto a = helical_params_.getStepRotationMatrix();
    vector3 r1(a[0], a[1], a[2]), r2(a[3], a[4], a[5]), r3(a[6], a[7], a[8]);
    step_rot = matrix3x3(r1, r2, r3);
    glbl_rot = helical_params_.getGlobalRotationMatrix();
    step_translate = helical_params_.getStepTranslationVec();
    glbl_translate = helical_params_.getGlobalTranslationVec();
    rng.seed(std::random_device()());
}

bool MonteCarloRotorSearch::run() {
    BaseUnit unit(base_a_, backbone_);
    Chain chain(unit, 10, runtime_params_.type);
    test_chain = chain.getChain();
    auto bu_a_mol = unit.getMol();
    auto bu_a_head_tail = unit.getBackboneLinkers();
    unsigned head = static_cast<unsigned>(bu_a_head_tail[0]),
             tail = static_cast<unsigned>(bu_a_head_tail[1]);


    bu_a_mol.Translate(glbl_translate);
    bu_a_mol.Rotate(glbl_rot.data());

    mol = OBMol(bu_a_mol);

    double* coords = bu_a_mol.GetCoordinates();
    uniform_real_distribution<double> dist(0, 2 * M_PI);
    uniform_real_distribution<double> one_zero_dist(0, 1);
    double k_effective = 0.59 / 5.15; // Angstroms^2
    monomer_num_coords = bu_a_mol.NumAtoms() * 3;


    OBRotorList rl;
    OBBitVec fix_bonds(backbone_.getMolecule().NumAtoms());
    auto base_indices = unit.getBaseIndexRange();
    for (unsigned i = static_cast<unsigned>(base_indices[0]); i <= base_indices[1]; ++i)
        fix_bonds.SetBitOn(i);

    rl.Setup(bu_a_mol);
    rl.SetFixAtoms(fix_bonds);
    rl.SetRotAtomsByFix(bu_a_mol);

    size_t search_size = runtime_params_.num_steps;
    OBRotorIterator ri;

    for (size_t search_index = 0; search_index < search_size; ++search_index) {
        OBRotor *r = rl.BeginRotor(ri);
        double best_dist = std::numeric_limits<double>::infinity(), cur_dist = best_dist;
        while (r) {

            bool accept = false;
            best_dist = std::numeric_limits<double>::infinity();
            while (!accept) {
                auto angle = dist(rng);
                r->SetToAngle(coords, angle);
                cur_dist = measureDistance(coords, head, tail);
                if (cur_dist < best_dist
                    || exp(-pow(cur_dist - best_dist, 2) / k_effective) > one_zero_dist(rng)) {
                    accept = true;
                    best_dist = cur_dist;
                }
            }
            r = rl.NextRotor(ri);
        }

//        double angle_stp = 2 * M_PI / search_size;
//        OBAtom *a = bu_a_mol.GetAtom(10), *b = bu_a_mol.GetAtom(5), *c = bu_a_mol.GetAtom(22), *d = bu_a_mol.GetAtom(18);
//        bu_a_mol.SetTorsion(a, b, c, d, angle_stp * search_index);
//        double cur_dist = measureDistance(coords, head, tail);
//        cout << "Index: " << search_index << ", distance: " << cur_dist << ", angle: " << angle_stp * search_index * RAD_TO_DEG << endl;

        // if accept, add to vector of coord_vec

        if (cur_dist < runtime_params_.max_distance) {

            auto data = chain.generateConformerData(coords, helical_params_);
//            cout << "Energy: " << data.total_energy;
//            cout << "others: " << data.totTorsionE << ", " << data.VDWE << ", " << data.torsionE << ", "
//                 << data.bondE << ", " << data.angleE << endl;
            if (!isPassingEFilter(data)) {
                delete data.coords;
            } else {
                data.monomer_coord = new double[monomer_num_coords];
                data.index = search_index;
                data.distance = cur_dist;
                memcpy(data.monomer_coord, coords, sizeof(double) * monomer_num_coords);
                print(data);
            }
        }
        if (search_index % 100000 == 0)
            cout << 100 * static_cast<double>(search_index) / search_size
                 << "%\t\tAccepted: " << conf_data_vec.size() << endl;
    }

    for (auto &v : conf_data_vec)
        delete v.monomer_coord;
    mol = bu_a_mol;

    return true;
}

double MonteCarloRotorSearch::measureDistance(double *coords, unsigned head, unsigned tail) {
    auto hi = 3 * (head - 1), ti = 3 * (tail - 1);
    vector3 head_coord(coords[hi], coords[hi + 1], coords[hi + 2]);
    vector3 tail_coord(coords[ti], coords[ti + 1], coords[ti + 2]);

    tail_coord *= step_rot;
    tail_coord += step_translate;
    return sqrt(head_coord.distSq(tail_coord));
}

bool MonteCarloRotorSearch::isPassingEFilter(const PNAB::ConformerData &conf_data) {
    vector<double> cur_vals = {conf_data.total_energy, conf_data.angleE, conf_data.bondE,
                               conf_data.VDWE, conf_data.totTorsionE};
    auto max_vals = runtime_params_.energy_filter;
    for (auto i = 0; i < max_vals.size(); ++i)
        if (max_vals[i] < cur_vals[i])
            return false;
    return true;
}

void MonteCarloRotorSearch::print(PNAB::ConformerData conf_data) {
    if (!conf_data.chain_coords_present) {
        cerr << "Trying to print conformer with no chain coordinates. Exiting..." << endl;
        exit(1);
    }
    ostringstream strs;
    filebuf fb;
    ofstream csv;

    // Set output format
    conv.SetOutFormat("PDB");

    // Clear string stream
    strs.str(std::string());

    // All conformers are named "conformer" + index
    strs << "conformer_" << conf_data.index << ".pdb";

    // Open file for writing...
    fb.open(strs.str().c_str(), std::ios::out);
    ostream fileStream(&fb);

    // Set the conformer save to file
    test_chain.SetCoordinates(conf_data.coords);
    conv.SetOutStream(&fileStream);
    conv.Write(&test_chain);
    fb.close();

    delete conf_data.coords;
    conf_data.chain_coords_present = false;
    conf_data_vec.push_back(conf_data);

    // csv file containing energy headings
    csv.open("energy_data.csv");
    csv << "Conformer Index, Energy (kcal/mol), Distance (A), Bond Energy, Angle Energy,";
    csv << " Torsion Energy, VDW Energy, Total Torsion Energy, RMSD (A)" << endl;
    std::sort(conf_data_vec.begin(), conf_data_vec.end());

    double *ref = conf_data_vec[0].monomer_coord;

    for (auto &v : conf_data_vec) {
        v.rmsd = calcRMSD(ref, v.monomer_coord, monomer_num_coords);
        csv << v.index  << ", " << v.total_energy << ", " << v.distance << ", " << v.bondE << ", "
            << v.angleE << ", " << v.torsionE     << ", " << v.VDWE     << ", " << v.totTorsionE << ", "
            << v.rmsd   << endl;
    }
    csv.close();
}
