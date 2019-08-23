//
// Created by jbarnett8 on 2/26/18.
//

#include <iomanip>
#include "ConformationSearch.h"

using namespace PNAB;
using namespace std;
using namespace OpenBabel;

ConformationSearch::ConformationSearch(RuntimeParameters &runtime_params, Backbone backbone,
                                       HelicalParameters &helical_params, Bases bases, string prefix) {
    prefix_ = prefix;
    runtime_params_ = runtime_params;
    strand_ = runtime_params_.strand;
    for (auto &v : strand_)
        transform(v.begin(), v.end(), v.begin(), ::tolower);
    auto name = strand_[0];
    transform(name.begin(), name.end(), name.begin(), ::tolower);
    base_a_ = bases.getBaseFromName(name);
    helical_params_ = helical_params;
    backbone_ = backbone;
    step_rot_ = helical_params_.getStepRotationOBMatrix(1);
    glbl_rot_ = helical_params_.getGlobalRotationOBMatrix();
    step_translate_ = helical_params_.getStepTranslationVec(1);
    glbl_translate_ = helical_params_.getGlobalTranslationVec();
    bases_ = bases;
    rng_.seed(std::random_device()());
    is_double_stranded_ = runtime_params_.is_double_stranded;
    is_hexad_ = runtime_params_.is_hexad;
    ff_type_ = runtime_params_.type;
}

std::string ConformationSearch::run() {
    // Choose search algorithm
    string search = runtime_params_.search_algorithm;
    transform(search.begin(), search.end(), search.begin(), ::tolower);

    if(search.find("weighted random search") != std::string::npos)
        return ConformationSearch::RandomSearch(true);

    else if (search.find("random search") != std::string::npos)
        return ConformationSearch::RandomSearch(false);

    else if (search.find("systematic search") != std::string::npos)
        return ConformationSearch::SystematicSearch();

    else if (search.find("weighted monte carlo search") != std::string::npos)
        return ConformationSearch::MonteCarloSearch(true);

    else if (search.find("monte carlo search") != std::string::npos)
        return ConformationSearch::MonteCarloSearch(false);

    else {
        cerr << search << " is unrecognized search algorithm" << endl;
        exit(1);
    }
}


std::string ConformationSearch::RandomSearch(bool weighted) {
    // Random search: weighted or uniform probability distribution

    // Get a base unit
    BaseUnit unit(base_a_, backbone_);
    auto range = unit.getBackboneIndexRange();
    backbone_range_ = {static_cast<unsigned >(range[0]), static_cast<unsigned >(range[1])};
    Chain chain(bases_, backbone_, strand_, ff_type_, backbone_range_, is_double_stranded_, is_hexad_, runtime_params_.strand_orientation);
    test_chain_ = chain.getChain();
    auto bu_a_mol = unit.getMol();
    auto bu_a_head_tail = unit.getBackboneLinkers();
    auto head = static_cast<unsigned>(bu_a_head_tail[0]),
         tail = static_cast<unsigned>(bu_a_head_tail[1]);
    auto fixed_bonds = unit.getFixedBonds();

    std::string output_string;

    double* coords = bu_a_mol.GetCoordinates();
    monomer_num_coords_ = bu_a_mol.NumAtoms() * 3;

    // Setup the rotor list
    OBRotorList rl;
    OBBitVec fix_bonds(backbone_.getMolecule().NumAtoms());
    auto base_indices = unit.getBaseIndexRange();
    for (unsigned i = static_cast<unsigned>(base_indices[0]); i <= base_indices[1]; ++i)
        fix_bonds.SetBitOn(i);
    rl.Setup(bu_a_mol);
    rl.SetFixAtoms(fix_bonds);
    rl.SetRotAtomsByFix(bu_a_mol);

    size_t search_size = runtime_params_.num_steps;

    // store all rotors in a vector
    OBRotorIterator ri;
    vector<OBRotor*> rotor_vector;
    OBRotor *r = rl.BeginRotor(ri);
    while(r) {
        bool save = true;
        unsigned a1 = r->GetDihedralAtoms()[1];
        unsigned a2 = r->GetDihedralAtoms()[2];
        for (auto f: fixed_bonds) {
            if ((a1 == f[0] && a2 == f[1]) || (a1 == f[1] && a2 == f[0])) {
                save = false;
            }
        }
        if (save)
            rotor_vector.push_back(r);
        r = rl.NextRotor(ri);
    }

    uniform_real_distribution<double> dist;
    vector <std::piecewise_linear_distribution<double>> dist_vector;

    // Get probability distribution
    if (weighted) {
        dist_vector = weighted_distributions(rotor_vector, bu_a_mol);
    }

    else {
        dist = uniform_real_distribution<double>(0, 2 * M_PI);
    }

    for (size_t search_index = 0; search_index < search_size; ++search_index) {

        // For random search, we rotate all dihedrals at every step 
        for (int i = 0; i < rotor_vector.size(); i++) {
            auto r = rotor_vector[i];
            double angle;
            // Choose a random angle
            if (weighted) {
                angle = dist_vector[i](rng_);
            }
            else {
                angle = dist(rng_);
            }
            r->SetToAngle(coords, angle);
        }

        // measure distance
        double cur_dist = measureDistance(coords, head, tail);

        // if accept, add to vector of coord_vec_
        if (cur_dist < runtime_params_.max_distance) {

            // Generate chain and compute energies; check whether energies are less than thresholds
            auto data = chain.generateConformerData(coords, helical_params_, runtime_params_.energy_filter);

            if (!data.accepted)
                delete[] data.coords;

            else {
                data.monomer_coord = new double[monomer_num_coords_];
                data.index = search_index;
                data.distance = cur_dist;
                memcpy(data.monomer_coord, coords, sizeof(double) * monomer_num_coords_);
                output_string = print(data);
            }
        }
        if (search_index % 100000 == 0) {
            printProgress(search_index, search_size);
        }
    }

    for (auto &v : conf_data_vec_)
        delete[] v.monomer_coord;

    return output_string;
}

std::string ConformationSearch::MonteCarloSearch(bool weighted) {
    // Monte Carlo search: weighted or uniform probability distribution

    // Get a base unit
    BaseUnit unit(base_a_, backbone_);
    auto range = unit.getBackboneIndexRange();
    backbone_range_ = {static_cast<unsigned >(range[0]), static_cast<unsigned >(range[1])};
    Chain chain(bases_, backbone_, strand_, ff_type_, backbone_range_, is_double_stranded_, is_hexad_, runtime_params_.strand_orientation);
    test_chain_ = chain.getChain();
    auto bu_a_mol = unit.getMol();
    auto bu_a_head_tail = unit.getBackboneLinkers();
    auto head = static_cast<unsigned>(bu_a_head_tail[0]),
         tail = static_cast<unsigned>(bu_a_head_tail[1]);
    auto fixed_bonds = unit.getFixedBonds();

    std::string output_string;

    double* coords = bu_a_mol.GetCoordinates();
    monomer_num_coords_ = bu_a_mol.NumAtoms() * 3;

    // Set up the rotor list
    OBRotorList rl;
    OBBitVec fix_bonds(backbone_.getMolecule().NumAtoms());
    auto base_indices = unit.getBaseIndexRange();
    for (unsigned i = static_cast<unsigned>(base_indices[0]); i <= base_indices[1]; ++i)
        fix_bonds.SetBitOn(i);
    rl.Setup(bu_a_mol);
    rl.SetFixAtoms(fix_bonds);
    rl.SetRotAtomsByFix(bu_a_mol);

    size_t search_size = runtime_params_.num_steps;

    // store all rotors in a vector
    OBRotorIterator ri;
    vector<OBRotor*> rotor_vector;
    OBRotor *r = rl.BeginRotor(ri);
    while(r) {
        bool save = true;
        unsigned a1 = r->GetDihedralAtoms()[1];
        unsigned a2 = r->GetDihedralAtoms()[2];
        for (auto f: fixed_bonds) {
            if ((a1 == f[0] && a2 == f[1]) || (a1 == f[1] && a2 == f[0])) {
                save = false;
            }
        }
        if (save)
            rotor_vector.push_back(r);
        r = rl.NextRotor(ri);
    }

    // Get probability distribution
    uniform_real_distribution<double> dist;
    vector <std::piecewise_linear_distribution<double>> dist_vector;

    if (weighted) {
        dist_vector = weighted_distributions(rotor_vector, bu_a_mol);
    }

    else {
        dist = uniform_real_distribution<double>(0, 2 * M_PI);
    }

    uniform_real_distribution<double> one_zero_dist(0, 1);
    double k_effective = 0.593 / 5.15; // Angstroms^2; 0.593 is kbT at 298 K in kcal/mol
    double best_dist = std::numeric_limits<double>::infinity();

    for (size_t search_index = 0; search_index < search_size; ++search_index) {
        // In Monte Carlo search, we randomly pick one dihedral and rotate it
        int index = rand()%rotor_vector.size();
        r = rotor_vector[index];

        // Save previous angle, in case the step is not good
        double old_angle = r->CalcTorsion(coords);

        // Choose a random angle
        double angle;
        if (weighted) {
            angle = dist_vector[index](rng_);
        }
        else {
            angle = dist(rng_);
        }
        r->SetToAngle(coords, angle);
        
        // measure distance
        double cur_dist = measureDistance(coords, head, tail);
        // Accept step if the new distance is less than the previous distance
        if (cur_dist < best_dist) {
            best_dist = cur_dist;
        }

        // If the new distance is larger than the previous distance,
        // accept step conditionally
        else if (one_zero_dist(rng_) < exp(-pow(cur_dist, 2) / k_effective)) {
            best_dist = cur_dist;
        }

        // If the step is not accepted, set the dihedral angle back to the previous state
        else {
            r->SetToAngle(coords, old_angle);
        }

        // if accept, add to vector of coord_vec_
        if (cur_dist < runtime_params_.max_distance) {

            // Generate chain and compute energies; check whether energies are less than thresholds
            auto data = chain.generateConformerData(coords, helical_params_, runtime_params_.energy_filter);

            if (!data.accepted)
                delete[] data.coords;

            else {
                data.monomer_coord = new double[monomer_num_coords_];
                data.index = search_index;
                data.distance = cur_dist;
                memcpy(data.monomer_coord, coords, sizeof(double) * monomer_num_coords_);
                output_string = print(data);
            }
        }
        if (search_index % 100000 == 0) {
            printProgress(search_index, search_size);
        }
    }

    for (auto &v : conf_data_vec_)
        delete[] v.monomer_coord;

    return output_string;
}


std::vector <std::piecewise_linear_distribution<double>> ConformationSearch::weighted_distributions(vector<OBRotor*> &rotor_vector, OBMol &bu_a_mol) {
    // Compute the weighted distributions bases on the energies of the dihedral angles

    // A vector to store the weighted distribution for each dihedral
    vector <std::piecewise_linear_distribution<double>> dist_vector;

    double kT = 0.593; // kbT at 298 K in kcal/mol

    // Set up rotor
    for (auto r: rotor_vector) {
        // For each dihedral angle, determine all atoms that are connected to the middle rotatable bond in the dihedral
        set<int> atoms;
        vector<int> v = r->GetDihedralAtoms();
        FOR_NBORS_OF_ATOM(nbr, bu_a_mol.GetAtom(v[1])) {
            if (nbr->GetIdx() == v[2])
                continue;
            FOR_NBORS_OF_ATOM(nbr2, bu_a_mol.GetAtom(v[2])) {
                if (nbr2->GetIdx() == v[1])
                    continue;
                atoms.insert(nbr->GetIdx());
                atoms.insert(v[1]);
                atoms.insert(v[2]);
                atoms.insert(nbr2->GetIdx());
        
            }
        }

        // Setup constraints to ignore all atoms except those that are connected to the rotatable bond
        OBFFConstraints constraints;
        FOR_ATOMS_OF_MOL(a, bu_a_mol) {
            if (!(atoms.find(a->GetIdx()) != atoms.end()))
                constraints.AddIgnore(a->GetIdx());
        }
        OBForceField *pFF = OBForceField::FindForceField(ff_type_);
        bool isKCAL_ = pFF->GetUnit().find("kcal") != string::npos;

        // Scan the energy of the dihedral angle at 1 degrees step size.
        double step = 1.0*M_PI/180.0, torsion = 0.0, sum = 0.0;
        vector<double> factors, intervals;

        for (int i=0; i < int(2*M_PI/step)+1; i++) {
            // save borders of intervals
            intervals.push_back(torsion);
            r->SetToAngle(bu_a_mol.GetCoordinates(), torsion);

            // Setup the force field
            pFF->Setup(bu_a_mol);
            pFF->SetConstraints(constraints);

            // Compute torsional energy
            double torsionE = pFF->E_Torsion(false);
            if (!isKCAL_)
                torsionE *= KJ_TO_KCAL;
            // Compute the Boltzmann factor at each angle and the sum
            double factor = exp(-1 * torsionE / kT);
            factors.push_back(factor);
            sum += factor;
            torsion += step;
        }


        // Compute the weights for each dihedral angle
        vector<double> weights;
        for (auto f: factors) {
            weights.push_back(f/sum);
        }
        
        // Generate the distribution. Intervals are weighted based on the energies.
        // Probabilities inside each interval are based on linear interpolation of the probabilities
        // at the border of the interval
        std::piecewise_linear_distribution<double> dist(intervals.begin(),intervals.end(),weights.begin());
        dist_vector.push_back(dist);
    }
    return dist_vector;
}


std::string ConformationSearch::SystematicSearch() {
    // Systematic search: Rotate each one of the dihedral angles at a given step size

    // Get a base unit
    BaseUnit unit(base_a_, backbone_);
    auto range = unit.getBackboneIndexRange();
    backbone_range_ = {static_cast<unsigned >(range[0]), static_cast<unsigned >(range[1])};
    Chain chain(bases_, backbone_, strand_, ff_type_, backbone_range_, is_double_stranded_, is_hexad_, runtime_params_.strand_orientation);
    test_chain_ = chain.getChain();
    auto bu_a_mol = unit.getMol();
    auto bu_a_head_tail = unit.getBackboneLinkers();
    auto head = static_cast<unsigned>(bu_a_head_tail[0]),
         tail = static_cast<unsigned>(bu_a_head_tail[1]);
    auto fixed_bonds = unit.getFixedBonds();

    std::string output_string;

    double* coords = bu_a_mol.GetCoordinates();
    monomer_num_coords_ = bu_a_mol.NumAtoms() * 3;

    // Setup the rotor list
    OBRotorList rl;
    OBBitVec fix_bonds(backbone_.getMolecule().NumAtoms());
    auto base_indices = unit.getBaseIndexRange();
    for (unsigned i = static_cast<unsigned>(base_indices[0]); i <= base_indices[1]; ++i)
        fix_bonds.SetBitOn(i);
    rl.Setup(bu_a_mol);
    rl.SetFixAtoms(fix_bonds);
    rl.SetRotAtomsByFix(bu_a_mol);

    // Determine the step size and the number of steps
    double dihedral_step = runtime_params_.dihedral_step*M_PI/180.0;
    size_t steps_per_bond = (size_t) round(2*M_PI/dihedral_step);
    size_t search_size = (size_t) pow(steps_per_bond, rl.Size());

    OBRotorIterator ri;

    OBRotor *r = rl.BeginRotor(ri);

    // Start by setting all angles to zero and set rotor vector
    vector<OBRotor*> rotor_vector;
    while(r) {
        bool save = true;
        unsigned a1 = r->GetDihedralAtoms()[1];
        unsigned a2 = r->GetDihedralAtoms()[2];
        for (auto f: fixed_bonds) {
            if ((a1 == f[0] && a2 == f[1]) || (a1 == f[1] && a2 == f[0])) {
                save = false;
            }
        }
        if (save) {
            rotor_vector.push_back(r);
            r->SetToAngle(coords, 0.0);
        }
        r = rl.NextRotor(ri);
    }

    for (size_t search_index = 1; search_index < search_size + 1; ++search_index) {
        // Choose rotor to change
        int rotor_index = 0;
        for (int i = 1; i < rotor_vector.size(); i++){
            size_t d = pow(steps_per_bond, i);
            int accept = search_index%d;
            if (accept==0) {
                rotor_index = i;
            }
        }

        // Set lower rotors to zero and start rotating
        for (int i=0; i < rotor_index; i++) {
            r = rotor_vector[i];
            r->SetToAngle(coords, 0.0);
        }

        // Set dihedral angle
        r = rotor_vector[rotor_index];
        double angle = r->CalcTorsion(coords) + dihedral_step;
        r->SetToAngle(coords, angle);

        double cur_dist = measureDistance(coords, head, tail);

        // if accept, add to vector of coord_vec_
        if (cur_dist < runtime_params_.max_distance) {

            auto data = chain.generateConformerData(coords, helical_params_, runtime_params_.energy_filter);

            if (!data.accepted)
                delete[] data.coords;

            else {
                data.monomer_coord = new double[monomer_num_coords_];
                data.index = search_index;
                data.distance = cur_dist;
                memcpy(data.monomer_coord, coords, sizeof(double) * monomer_num_coords_);
                output_string = print(data);
            }
        }
        if (search_index % 100000 == 0) {
            printProgress(search_index, search_size);
        }
    }

    for (auto &v : conf_data_vec_)
        delete[] v.monomer_coord;

    return output_string;
}

void ConformationSearch::printProgress(std::size_t search_index, std::size_t search_size) {
    // Print progress
    auto prgrs = conf_data_vec_.size();
    cout << setw(8);
    cout << 100 * static_cast<double>(search_index) / search_size;
    cout << "%\tAccepted: " << setw(8) << prgrs;
    if (prgrs > 0) {
        cout << ", Best Conformer (distance, energy): (" << setw(10) << conf_data_vec_[0].distance
             << ", " << setw(10) << conf_data_vec_[0].total_energy << ") -- " << prefix_ << "_"
             << conf_data_vec_[0].index << ".pdb" << endl;
    } else {
        cout << endl;
    }

}

double ConformationSearch::measureDistance(double *coords, unsigned head, unsigned tail) {
    // Measure distance between atom linkers in two consecutive backbones
    auto hi = 3 * (head - 1), ti = 3 * (tail - 1);
    vector3 head_coord(coords[hi], coords[hi + 1], coords[hi + 2]);
    vector3 tail_coord(coords[ti], coords[ti + 1], coords[ti + 2]);

    tail_coord += glbl_translate_;
    tail_coord *= glbl_rot_;
    head_coord += glbl_translate_;
    head_coord *= glbl_rot_;
    tail_coord += step_translate_;
    tail_coord *= step_rot_;
    return sqrt(head_coord.distSq(tail_coord));
}

std::string ConformationSearch::print(PNAB::ConformerData conf_data) {
    // Save data

    if (!conf_data.chain_coords_present) {
        cerr << "Trying to print conformer with no chain_ coordinates. Exiting..." << endl;
        exit(1);
    }
    ostringstream strs;
    filebuf fb;
    stringstream output_string;

    // Set output format
    conv_.SetOutFormat("PDB");

    // Clear string stream
    strs.str(std::string());

    // All conformers are named prefix_ + "_" + index
    strs << prefix_ << "_" << conf_data.index << ".pdb";

    // Open file for writing...
    fb.open(strs.str().c_str(), std::ios::out);
    ostream fileStream(&fb);

    // Set the conformer save to file
    test_chain_.SetTitle(strs.str().c_str());
    test_chain_.SetCoordinates(conf_data.coords);
    conv_.SetOutStream(&fileStream);
    conv_.Write(&test_chain_);
    fb.close();

    delete[] conf_data.coords;
    conf_data.chain_coords_present = false;
    conf_data_vec_.push_back(conf_data);

    output_string << "# Prefix, Conformer Index, Distance (A), Bond Energy, Angle Energy, Torsion Energy, VDW Energy, Total Energy, Fixed Torsion Energy (kcal/mol)";
    output_string << "RMSD" << endl;
    std::sort(conf_data_vec_.begin(), conf_data_vec_.end());

    double *ref = conf_data_vec_[0].monomer_coord;

    for (auto &v : conf_data_vec_) {
        v.rmsd = calcRMSD(ref, v.monomer_coord, monomer_num_coords_);
        output_string << prefix_ << ", " << v.index  << ", " << v.distance << ", " << v.bondE << ", " << v.angleE << ", "
            << v.torsionE << ", " << v.VDWE << ", " << v.total_energy << ", " << v.fixed_torsionE << ", " << v.rmsd << endl;
    }


    return output_string.str();
}
