/**@file
 * @brief A file for defining various search algorithm functions
 */

#include <iomanip>
#include <set>
#include "ConformationSearch.h"

using namespace PNAB;
using namespace std;
using namespace OpenBabel;

ConformationSearch::ConformationSearch(RuntimeParameters &runtime_params, Backbone &backbone,
                                       HelicalParameters &helical_params, Bases bases, string prefix, bool verbose) {

    // Set parameters
    runtime_params_ = runtime_params;
    backbone_ = backbone;
    helical_params_ = helical_params;
    bases_ = bases;
    prefix_ = prefix;
    verbose_ = verbose;

    // Update the names of the bases in the strand
    for (auto &v : runtime_params_.strand)
        transform(v.begin(), v.end(), v.begin(), ::tolower);
    // Get the first base in the strand for the test nucleotide
    auto name = runtime_params_.strand[0];
    transform(name.begin(), name.end(), name.begin(), ::tolower);
    auto base_a = bases.getBaseFromName(name);

    // Random number seed
    rng_.seed(runtime_params.seed);

    // Get the required translations and rotations
    step_rot_ = helical_params_.getStepRotationOBMatrix(1);
    glbl_rot_ = helical_params_.getGlobalRotationOBMatrix();
    step_translate_ = helical_params_.getStepTranslationVec(1);
    glbl_translate_ = helical_params_.getGlobalTranslationVec();

    // Get a base unit
    unit = BaseUnit(base_a, backbone_);
    auto range = unit.getBackboneIndexRange();
    backbone_range_ = {static_cast<unsigned >(range[0]), static_cast<unsigned >(range[1])};
    bu_a_mol = unit.getMol();
    auto bu_a_head_tail = unit.getBackboneLinkers();
    head = static_cast<unsigned>(bu_a_head_tail[0]);
    tail = static_cast<unsigned>(bu_a_head_tail[1]);
    auto fixed_bonds = unit.getFixedBonds();
    coords = bu_a_mol.GetCoordinates();
    monomer_num_coords_ = bu_a_mol.NumAtoms() * 3;

    // Setup the rotor list
    OBBitVec fix_bonds(backbone_.getMolecule().NumAtoms());
    // fix the base atoms
    auto base_indices = unit.getBaseIndexRange();
    for (unsigned i = static_cast<unsigned>(base_indices[0]); i <= base_indices[1]; ++i)
        fix_bonds.SetBitOn(i);

    rl.Setup(bu_a_mol);
    rl.SetFixAtoms(fix_bonds);
    rl.SetRotAtomsByFix(bu_a_mol);

    // store all rotors in a vector
    OBRotorIterator ri;
    OBRotor *r = rl.BeginRotor(ri);
    while(r) {
        bool save = true;
        unsigned a1 = r->GetDihedralAtoms()[1];
        unsigned a2 = r->GetDihedralAtoms()[2];
        // If the dihedral angle is fixed, do not add it to the rotor list
        for (auto f: fixed_bonds) {
            if ((a1 == f[0] && a2 == f[1]) || (a1 == f[1] && a2 == f[0])) {
                save = false;
            }
        }
        if (save)
            rotor_vector.push_back(r);
        r = rl.NextRotor(ri);
    }

    // Header
    std::string header = "# Prefix, Conformer Index, Distance (Angstroms), Bond Energy (kcal/mol), Angle Energy (kcal/mol), "
                         "Torsion Energy (kcal/mol/nucleotide), Van der Waals Energy (kcal/mol/nucleotide), "
                         "Total Energy (kcal/mol/nucleotide)";
    for (int i=0; i < rotor_vector.size(); i++) {
        header += ", Dihedral " + to_string(i+1) + " (degrees)";
    }

    output_string += header + "\n";
}

std::string ConformationSearch::run() {

    // Choose search algorithm
    string search = runtime_params_.search_algorithm;
    transform(search.begin(), search.end(), search.begin(), ::tolower);

    if(search.find("weighted random search") != std::string::npos)
        ConformationSearch::RandomSearch(true);

    else if (search.find("random search") != std::string::npos)
        ConformationSearch::RandomSearch(false);

    else if (search.find("systematic search") != std::string::npos)
        ConformationSearch::SystematicSearch();

    else if (search.find("weighted monte carlo search") != std::string::npos)
        ConformationSearch::MonteCarloSearch(true);

    else if (search.find("monte carlo search") != std::string::npos)
        ConformationSearch::MonteCarloSearch(false);

    else if (search.find("genetic algorithm search") != std::string::npos)
        ConformationSearch::GeneticAlgorithmSearch();

    else {
        throw std::runtime_error(search + " is unrecognized search algorithm");
    }

    // return the CSV output string
    return output_string;
}


void ConformationSearch::GeneticAlgorithmSearch() {
    // Genetic algorithm search

    // Setup chain
    Chain chain(bases_, backbone_, runtime_params_.strand, runtime_params_.ff_type, backbone_range_,
                runtime_params_.is_hexad, runtime_params_.build_strand, runtime_params_.strand_orientation);

    // Set the search size; the number of generations in the genetic algorithm search
    size_t search_size = runtime_params_.num_steps;

    // Use uniform probability between 0 and 2 pi for the dihedral angles 
    uniform_real_distribution<double> dist = uniform_real_distribution<double>(0, 2 * M_PI);
    // Use uniform probability between 0 and 1 for the use in the crossover and mutation 
    uniform_real_distribution<double> dist2 = uniform_real_distribution<double>(0, 1);
    uniform_int_distribution<int> selection(0, rotor_vector.size() - 1);

    // Genetic algorithm parameters
    int numConformers = runtime_params_.population_size; // Population size
    int elites = 0; // Number of intact survivors; Not used now. Set to zero 
    double mutation_rate = runtime_params_.mutation_rate; // Mutation rate; large values prevent from getting stuck
    double crossover_rate = runtime_params_.crossover_rate; // Mating rate
    vector<pair<double, vector<double>>> population; //(fitness (1/distance), vector of the torsional angles)

    // Initialize population with random dihedral angles
    for (int con=0; con < numConformers; con++ ) {
        vector<double> state;
        for (int i=0; i < rotor_vector.size(); i++) {
            // Set all the dihedral angles to random values
            state.push_back(dist(rng_));
            auto r = rotor_vector[i];
            r->SetToAngle(coords, state[i]);
        }
        // Compute the fitness
        double cur_dist = measureDistance(coords, head, tail);
        population.push_back(make_pair(1.0/cur_dist, state));
    }

    size_t save_index = 0; // Different index used to save pdb files
    vector<int> save_cur_dist; // Save distance to avoid processing identical conformers

    // Loop over the number of generations
    for (size_t search_index=0; search_index < search_size; search_index++) {

        // print progress roughly every 10%
        if (fmod(search_index, search_size/10) == 0 && verbose_) {
            printProgress(search_index, search_size);
        }

        // Sort population by fitness
        sort(population.begin(), population.end());
        vector<double> weights;
        for (int i=0; i < population.size(); i++) {
            // It seems rank-based weights are better in our case
            // Weighting based on the fitness value seems to be too strong
            weights.push_back(i + 1); // The weight is proportional to the rank
        }

        // Selection probability for choosing an individuals for crossover and mutation
        // Generate the index of the selected individual
        discrete_distribution<> selection_probability (weights.begin(), weights.end()); 

        // Vector to store new generation
        vector<pair<double, vector<double>>> new_generation;

        // Save some candidates for the next generation
        // Not implemented now
        //for(int i=0; i < elites; i++) {
        //    int r = selection_probability(rng_);
        //    new_generation.push_back(population[r]); 
        //}

        // Generate offsprings
        for (int i=0; i < (numConformers-elites)/2; i++) { // No survivors now

            // Select parents; allow self-mating
            int r1 = selection_probability(rng_);
            int r2 = selection_probability(rng_);
            vector<double> parent1 = population[r1].second;
            vector<double> parent2 = population[r2].second;

            // Set the offsprings to be equal to the parents and then perform mating and mutation
            vector<vector<double>> offspring = {parent1, parent2};

            // Crossover
            // Get a random number
            double v = dist2(rng_);
            if (v < crossover_rate) {
                // Exchange one dihedral angle for each offspring
                // Maybe a different choice for the mating procedure is preferred
                int index = selection(rng_)%offspring[0].size();
                offspring[0][index] = parent2[index];
                offspring[1][index] = parent1[index];
            }

            // Mutate
            // Get a random number
            v = dist2(rng_);
            if (v < mutation_rate) {
                // Select a random dihedral angle and change it
                // Maybe a different choice for the mutation procedure is preferred 
                int index1 = selection(rng_)%offspring[0].size();
                int index2 = selection(rng_)%offspring[1].size();
                offspring[0][index1] = dist(rng_);
                offspring[1][index2] = dist(rng_);
            }

            // Compute fitness (inverse distance) and see if the two offsprings are good candidates
            for (int kid=0; kid< 2; kid++) {
                save_index++;
                // Set the coordinates to the new angles
                for (int j=0; j < offspring[kid].size(); j++) {
                    auto r = rotor_vector[j];
                    r->SetToAngle(coords, offspring[kid][j]);
                }

                // Compute fitness
                double cur_dist = measureDistance(coords, head, tail);
                new_generation.push_back(make_pair(1.0/cur_dist, offspring[kid]));

                // if accept, add to vector of coord_vec_
                if (cur_dist < runtime_params_.max_distance) {
                    // Do not process previously explored structures
                    // We cannot prevent getting identical structures. But we do not need to compute energies for it
                    // or save it again
                    if (find(save_cur_dist.begin(), save_cur_dist.end(), int(cur_dist*1e6)) != save_cur_dist.end())
                        continue;
                    save_cur_dist.push_back(int(cur_dist*1e6));

                    // Generate chain and compute energies; check whether energies are less than thresholds
                    auto data = chain.generateConformerData(coords, helical_params_, runtime_params_.energy_filter);

                    if (data.accepted) {
                        // Save the candidate
                        data.monomer_coord = new double[monomer_num_coords_];
                        data.index = save_index;
                        data.distance = cur_dist;
                        memcpy(data.monomer_coord, coords, sizeof(double) * monomer_num_coords_);
                        reportData(data);
                        delete[] data.monomer_coord;
                    }
                }
            }
        }

        // Set the population to the new generation
        population = new_generation;

    }

    return;

}

void ConformationSearch::RandomSearch(bool weighted) {
    // Random search: weighted or uniform probability distribution

    // Setup chain
    Chain chain(bases_, backbone_, runtime_params_.strand, runtime_params_.ff_type, backbone_range_,
                runtime_params_.is_hexad, runtime_params_.build_strand, runtime_params_.strand_orientation);

    // Set the search size;
    size_t search_size = runtime_params_.num_steps;

    // Inititalize two different probability distibutions
    // We will choose based on whether we want weighted or unweighted search
    uniform_real_distribution<double> dist;
    vector <std::piecewise_linear_distribution<double>> dist_vector;

    // Get probability distribution
    if (weighted) {
        // Get the weighted distibution for each dihedral angle
        dist_vector = WeightedDistributions();
    }

    else {
        // Use a uniform distribution
        dist = uniform_real_distribution<double>(0, 2 * M_PI);
    }

    // Loop over the number of iterations
    for (size_t search_index = 0; search_index < search_size; ++search_index) {

        // print progress roughly every 10%
        if (fmod(search_index, search_size/10) == 0 && verbose_) {
            printProgress(search_index, search_size);
        }

        // For random search, we rotate all dihedrals at every step 
        for (int i = 0; i < rotor_vector.size(); i++) {
            // pick the rotor
            auto r = rotor_vector[i];
            double angle;
            // Choose a random angle
            if (weighted) {
                angle = dist_vector[i](rng_);
            }
            else {
                angle = dist(rng_);
            }
            // Set new angle
            r->SetToAngle(coords, angle);
        }

        // measure distance
        double cur_dist = measureDistance(coords, head, tail);

        // if accept, add to vector of coord_vec_
        if (cur_dist < runtime_params_.max_distance) {

            // Generate chain and compute energies; check whether energies are less than thresholds
            auto data = chain.generateConformerData(coords, helical_params_, runtime_params_.energy_filter);

            if (data.accepted) {
                // Save the candidate
                data.monomer_coord = new double[monomer_num_coords_];
                data.index = search_index;
                data.distance = cur_dist;
                memcpy(data.monomer_coord, coords, sizeof(double) * monomer_num_coords_);
                reportData(data);
                delete[] data.monomer_coord;
            }
        }
    }

    return;
}

void ConformationSearch::MonteCarloSearch(bool weighted) {
    // Monte Carlo search: weighted or uniform probability distribution

    // Setup chain
    Chain chain(bases_, backbone_, runtime_params_.strand, runtime_params_.ff_type, backbone_range_,
                runtime_params_.is_hexad, runtime_params_.build_strand, runtime_params_.strand_orientation);

    // Set the search size;
    size_t search_size = runtime_params_.num_steps;

    // Inititalize two different probability distibutions
    // We will choose based on whether we want weighted or unweighted search
    uniform_real_distribution<double> dist;
    vector <std::piecewise_linear_distribution<double>> dist_vector;

    // Get probability distribution
    if (weighted) {
        // Get the weighted distibution for each dihedral angle
        dist_vector = WeightedDistributions();
    }

    else {
        // Use a uniform distribution
        dist = uniform_real_distribution<double>(0, 2 * M_PI);
    }

    // Set a uniform random distribution for the acceptance and rejection step
    uniform_real_distribution<double> one_zero_dist(0, 1);
    

    double k = 300; // kcal/mol/Angstroms^2; Bond stretching force constant
                    // P-O5' bond stretching contsant in CHARMM36 for ON2-P bond: 270 kcal/mol/Angstrom^2
                    // V(bond) = Kb(b-b0)^2
                    // This is a rough value to get us close to the accepted distance
    double kT = runtime_params_.monte_carlo_temperature * BOLTZMANN; // kbT in kcal/mol

    // Set the initial distance to infinity
    double best_dist = std::numeric_limits<double>::infinity();

    // Loop over the number of iterations
    for (size_t search_index = 0; search_index < search_size; ++search_index) {

        // print progress roughly every 10%
        if (fmod(search_index, search_size/10) == 0 && verbose_) {
            printProgress(search_index, search_size);
        }

        // In Monte Carlo search, we randomly pick two dihedrals and rotate them
        int n_rotations = 2;
        // Store the indices of the rotors
        vector<int> indices;
        for (int i=0; i < rotor_vector.size(); i++)
            indices.push_back(i); 

        // Randomly choose two indices
        shuffle(indices.begin(), indices.end(), rng_);
        vector<int> rotated_indices = {indices[0], indices[1]};

        // Save the old angles and the rotated indices
        vector<double> old_angles;
        for (int i=0; i < n_rotations; i++) {
            // Pick the index
            int index = rotated_indices[i];
            
            // Get the rotor
            auto r = rotor_vector[index];

            // Save previous angle, in case the step is not good
            double old_angle = r->CalcTorsion(coords);
            old_angles.push_back(old_angle);

            // Choose a random angle
            double angle;
            if (weighted) {
                angle = dist_vector[index](rng_);
            }
            else {
                angle = dist(rng_);
            }
            // Set new angle
            r->SetToAngle(coords, angle);
        }
        // measure distance
        double cur_dist = measureDistance(coords, head, tail);
        // Accept step if the new distance is less than the previous distance
        if (cur_dist < best_dist) {
            best_dist = cur_dist;
        }

        // If the new distance is larger than the previous distance,
        // accept step conditionally
        else if (one_zero_dist(rng_) < exp(-k*pow((cur_dist-best_dist), 2) / kT)) {
            best_dist = cur_dist;
        }

        // If the step is not accepted, set the dihedral angle back to the previous state
        else {
            for (int i=0; i < n_rotations; i++) {
                auto r = rotor_vector[rotated_indices[i]];
                r->SetToAngle(coords, old_angles[i]);
            }
            // Go to the next iteration
            continue;
        }

        // if accept, add to vector of coord_vec_
        if (cur_dist < runtime_params_.max_distance) {

            // Generate chain and compute energies; check whether energies are less than thresholds
            auto data = chain.generateConformerData(coords, helical_params_, runtime_params_.energy_filter);

            if (data.accepted) {
                // Save the candidate
                data.monomer_coord = new double[monomer_num_coords_];
                data.index = search_index;
                data.distance = cur_dist;
                memcpy(data.monomer_coord, coords, sizeof(double) * monomer_num_coords_);
                reportData(data);
                delete[] data.monomer_coord;
            }
        }
    }

    return;
}


std::vector <std::piecewise_linear_distribution<double>> ConformationSearch::WeightedDistributions() {
    // Compute the weighted distributions bases on the energies of the dihedral angles

    // A vector to store the weighted distribution for each dihedral
    vector <std::piecewise_linear_distribution<double>> dist_vector;

    // The temperatue used in weighting the dihedral angles
    double kT = runtime_params_.weighting_temperature * BOLTZMANN; // kbT in kcal/mol

    // Set up rotor
    for (auto r: rotor_vector) {
        // For each dihedral angle, determine all atoms that are connected to the rotatable bond in the dihedral
        vector<int> v = r->GetDihedralAtoms();
        set<int> atoms = {v[1], v[2]};
        // Add neighbors of the first atom
        FOR_NBORS_OF_ATOM(nbr, bu_a_mol.GetAtom(v[1])) {
            if (nbr->GetIdx() == v[2])
                continue;
            atoms.insert(nbr->GetIdx());
        }
        // Add neighbors of the second atom
        FOR_NBORS_OF_ATOM(nbr2, bu_a_mol.GetAtom(v[2])) {
            if (nbr2->GetIdx() == v[1])
                continue;
            atoms.insert(nbr2->GetIdx());
        }

        // Setup constraints to ignore all atoms except those that are connected to the rotatable bond
        OBFFConstraints constraints;
        FOR_ATOMS_OF_MOL(a, bu_a_mol) {
            if (!(atoms.find(a->GetIdx()) != atoms.end()))
                constraints.AddIgnore(a->GetIdx());
        }
        // Setup force field
        OBForceField *pFF = OBForceField::FindForceField(runtime_params_.ff_type);
        if (!pFF) {
            throw std::runtime_error("Cannot find force field.");
        }
        // Check whether the energies are computed in kcal/mol
        bool isKCAL_ = pFF->GetUnit().find("kcal") != string::npos;

        // Scan the energy of the dihedral angle at 5 degrees step size.
        double step = 5.0*M_PI/180.0, torsion = 0.0, sum = 0.0;
        // The weighting factor and the intervals for each dihedral
        vector<double> factors, intervals;

        for (int i=0; i < int(2*M_PI/step)+1; i++) {
            // save borders of intervals
            intervals.push_back(torsion);
            // Set to the torsion
            r->SetToAngle(bu_a_mol.GetCoordinates(), torsion);

            // Setup the force field
            pFF->Setup(bu_a_mol);
            pFF->SetConstraints(constraints);

            // Compute torsional energy
            double torsionE = pFF->E_Torsion(false);
            if (!isKCAL_)
                torsionE *= KJ_TO_KCAL;
            // Compute the Boltzmann factor at each angle and the Boltzmann sum
            double factor = exp(-1 * torsionE / kT);
            factors.push_back(factor);
            sum += factor;
            // Increase the torsion
            torsion += step;
        }


        // Compute the weights for each dihedral angle
        // Divide the Boltzmann factor of each dihedral angle by the Boltzmann sum
        vector<double> weights;
        for (auto f: factors) {
            weights.push_back(f/sum);
        }
        
        // Generate the distribution. Intervals are weighted based on the Boltzmann ratio.
        // Probabilities inside each interval are based on linear interpolation of the probabilities
        // at the border of the interval
        std::piecewise_linear_distribution<double> dist(intervals.begin(),intervals.end(),weights.begin());
        dist_vector.push_back(dist);
    }

    return dist_vector;
}


void ConformationSearch::SystematicSearch() {
    // Systematic search: Rotate each one of the dihedral angles at a given step size

    // Setup chain
    Chain chain(bases_, backbone_, runtime_params_.strand, runtime_params_.ff_type, backbone_range_,
                runtime_params_.is_hexad, runtime_params_.build_strand, runtime_params_.strand_orientation);

    // Determine the step size and the number of steps
    // The number of steps is (360/dihedral_step)^(number of rotors)
    double dihedral_step = runtime_params_.dihedral_step*M_PI/180.0;
    size_t steps_per_bond = (size_t) round(2*M_PI/dihedral_step);
    size_t search_size = (size_t) pow(steps_per_bond, rotor_vector.size());

    // Set all rotors to zero to start
    for (auto r: rotor_vector)
        r->SetToAngle(coords, 0.0);

    // Loop over the number of iterations
    for (size_t search_index = 1; search_index < search_size + 1; ++search_index) {

        // print progress roughly every 10%
        if (fmod(search_index, search_size/10) == 0 && verbose_) {
            printProgress(search_index, search_size);
        }

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
            auto r = rotor_vector[i];
            r->SetToAngle(coords, 0.0);
        }

        // Set dihedral angle
        auto r = rotor_vector[rotor_index];
        double angle = r->CalcTorsion(coords) + dihedral_step;
        r->SetToAngle(coords, angle);

        // measure distance
        double cur_dist = measureDistance(coords, head, tail);

        // if accept, add to vector of coord_vec_
        if (cur_dist < runtime_params_.max_distance) {
            // Generate chain and compute energies; check whether energies are less than thresholds
            auto data = chain.generateConformerData(coords, helical_params_, runtime_params_.energy_filter);

            if (data.accepted) {
                // Save the candidate
                data.monomer_coord = new double[monomer_num_coords_];
                data.index = search_index;
                data.distance = cur_dist;
                memcpy(data.monomer_coord, coords, sizeof(double) * monomer_num_coords_);
                reportData(data);
                delete[] data.monomer_coord;
            }
        }
    }

    return;
}

void ConformationSearch::printProgress(std::size_t search_index, std::size_t search_size) {
    // Print progress
    cout << prefix_ << ": ";
    cout << 100 * static_cast<double>(search_index) / search_size;
    cout << "%\tAccepted: " << number_of_candidates << endl;
}

double ConformationSearch::measureDistance(double *coords, unsigned head, unsigned tail) {
    // Measure distance between atom linkers in two consecutive backbones
    auto hi = 3 * (head - 1), ti = 3 * (tail - 1);
    vector3 head_coord(coords[hi], coords[hi + 1], coords[hi + 2]);
    vector3 tail_coord(coords[ti], coords[ti + 1], coords[ti + 2]);

    // Perform the global rotation and translation for both the head and the tail
    tail_coord += glbl_translate_;
    tail_coord *= glbl_rot_;
    head_coord += glbl_translate_;
    head_coord *= glbl_rot_;
    // Perform the step translation and rotation to get the coordinate for the tail
    // atom in the next nucleotide
    tail_coord += step_translate_;
    tail_coord *= step_rot_;

    // Return the distance
    return sqrt(head_coord.distSq(tail_coord));
}

void ConformationSearch::reportData(PNAB::ConformerData &conf_data) {
    
    // Increase number of accepted candidates
    number_of_candidates += 1;

    // Setup variables for saving the structure of the accepted candidate
    ostringstream strs;
    filebuf fb;
    stringstream output_stringstream;

    // Set output format
    conv_.SetOutFormat("PDB");

    // Clear string stream
    strs.str(std::string());

    // All conformers are named prefix_ + "_" + index
    strs << prefix_ << "_" << conf_data.index << ".pdb";

    // Open file for writing...
    fb.open(strs.str().c_str(), std::ios::out);
    ostream fileStream(&fb);

    // Set conformer data and save to file
    conf_data.molecule.SetTitle(strs.str().c_str());

    OBPairData *pairdata = new OBPairData;
    pairdata->SetAttribute("AUTHOR");
    pairdata->SetValue("    The proto-Nucleic Acid Builder");
    conf_data.molecule.CloneData(pairdata);
    delete pairdata;

    vector<string> labels = {"Helical Rise (Angstroms)", "X-Displacement (Angstroms)", "Y-Displacement (Angstrom)",
                             "Helical Twist (degrees)", "Inclination (degrees)", "Tip (degrees)",
                             "Distance (Angstroms)", "Bond Energy (kcal/mol)", "Angle Energy (kcal/mol)", "Torsion Energy (kcal/mol/nucleotide)",
                             "Van der Waals Energy (kcal/mol/nucleotide)", "Total Energy (kcal/mol/nucleotide)"};
    vector<double> data = {helical_params_.h_rise, helical_params_.x_displacement, helical_params_.y_displacement,
                           helical_params_.h_twist, helical_params_.inclination, helical_params_.tip,
                           conf_data.distance, conf_data.bondE, conf_data.angleE, conf_data.torsionE, conf_data.VDWE, conf_data.total_energy};

    for (int i=0; i < rotor_vector.size(); i++) {
        labels.push_back("Dihedral " + to_string(i+1) + " (degrees)");
        double angle = rotor_vector[i]->CalcTorsion(conf_data.monomer_coord) * 180.0/M_PI;
        conf_data.dihedral_angles.push_back(angle);
        data.push_back(angle);
    }

    for (int i=0; i < 12 + rotor_vector.size(); i++) {
        OBPairData *pairdata = new OBPairData;
        pairdata->SetAttribute("TITLE");
        pairdata->SetValue("    " + labels[i] + ": " + to_string(data[i]));
        conf_data.molecule.CloneData(pairdata);
        delete pairdata;
    }

    conv_.SetOutStream(&fileStream);
    conv_.Write(&conf_data.molecule);
    fb.close();

    // Now we store the properties of the accepted candidates
    auto v = conf_data;
    output_stringstream << prefix_ << ", " << v.index  << ", " << v.distance << ", " << v.bondE << ", " << v.angleE << ", "
        << v.torsionE << ", " << v.VDWE << ", " << v.total_energy;
    for (auto angle: v.dihedral_angles)
        output_stringstream << ", " << angle;
    output_stringstream << endl;

    // Update the output string
    output_string += output_stringstream.str();

    return;
}
