/**@file
 * @brief A file for defining various methods for defining options
 */

#include "Containers.h"
#include <openbabel/math/align.h>

using namespace std;
using namespace PNAB;
using namespace OpenBabel;

void HelicalParameters::computeHelicalParameters() {
    if (!is_helical) {
        vector<vector3> ref_frame = StepParametersToReferenceFrame();
        ReferenceFrameToHelicalParameters(ref_frame[0], ref_frame[1], ref_frame[2], ref_frame[3]);
    }

    return;
}

vector3 HelicalParameters::getGlobalTranslationVec(bool is_base_pair, bool is_second_strand) {
    if (!is_base_pair)
        return vector3(x_displacement, y_displacement, 0);
    else if (!is_second_strand)
        return vector3(0.5*shear, 0.5*stretch, 0);
    else
        return vector3(-0.5*shear, -0.5*stretch, 0);
}

vector3 HelicalParameters::getStepTranslationVec(unsigned n, bool is_base_pair, bool is_second_strand) {
    if (!is_base_pair)
        return vector3(0, 0, n * h_rise);
    else if (!is_second_strand)
        return vector3(0, 0, 0.5*stagger);
    else
        return vector3(0, 0, -0.5*stagger);
}


matrix3x3 HelicalParameters::getGlobalRotationMatrix(bool is_base_pair, bool is_second_strand) {
    double angle1, angle2;
    if (!is_base_pair)
        angle1 = inclination * DEG_TO_RAD, angle2 = tip * DEG_TO_RAD;
    else if (!is_second_strand)
        angle1 = 0.5 * buckle * DEG_TO_RAD, angle2 = 0.5 * propeller * DEG_TO_RAD;
    else
        angle1 = -0.5 * buckle * DEG_TO_RAD, angle2 = -0.5 * propeller * DEG_TO_RAD;

    double Lambda = sqrt(angle1*angle1 + angle2*angle2);
    
    vector3 axis;
    if (Lambda != 0)
        axis = vector3(angle1/Lambda, angle2/Lambda, 0);
    else
       axis = vector3(0, 1, 0);

    matrix3x3 result = rodrigues_formula(axis, Lambda);

    return result;
}

matrix3x3 HelicalParameters::getStepRotationMatrix(unsigned n, bool is_base_pair, bool is_second_strand) {
    double angle;
    if (!is_base_pair)
        angle = n * h_twist * DEG_TO_RAD;
    else if (!is_second_strand)
        angle = 0.5 * opening * DEG_TO_RAD, n = 1;
    else
        angle = -0.5 * opening * DEG_TO_RAD, n = 1;

    matrix3x3 r_mat = matrix3x3(vector3(cos(angle), -sin(angle), 0),
                                vector3(sin(angle), cos(angle), 0),
                                vector3(0, 0, 1));

    return r_mat;
};

matrix3x3 HelicalParameters::rodrigues_formula(vector3 axis_vector, double theta) {

    array<double, 3> axis = {axis_vector.GetX(), axis_vector.GetY(), axis_vector.GetZ()};
    array<double, 9> m{};
    m[0] = cos(theta) + axis[0]*axis[0]*(1 - cos(theta));
    m[1] = axis[0]*axis[1]*(1 - cos(theta)) - axis[2]*sin(theta);
    m[2] = axis[1]*sin(theta) + axis[0]*axis[2]*(1 - cos(theta));
    m[3] = axis[2]*sin(theta) + axis[0]*axis[1]*(1 - cos(theta));
    m[4] = cos(theta) + axis[1]*axis[1]*(1 - cos(theta));
    m[5] = -axis[0]*sin(theta) + axis[1]*axis[2]*(1 - cos(theta));
    m[6] = -axis[1]*sin(theta) + axis[0]*axis[2]*(1 - cos(theta));
    m[7] = axis[0]*sin(theta) + axis[1]*axis[2]*(1 - cos(theta));
    m[8] = cos(theta) + axis[2]*axis[2]*(1 - cos(theta));

    matrix3x3 result = matrix3x3(vector3(m[0], m[1], m[2]),
                                 vector3(m[3], m[4], m[5]),
                                 vector3(m[6], m[7], m[8]));
    return result;
};


vector<vector3> HelicalParameters::StepParametersToReferenceFrame() {
    vector3 old_origin(0.0, 0.0, 0.0);
    vector3 x(1.0, 0.0, 0.0);
    vector3 y(0.0, 1.0, 0.0);
    vector3 z(0.0, 0.0, 1.0);

    // Compute the RollTill axis and angle; Equation 1
    vector3 hinge(tilt*DEG_TO_RAD, roll*DEG_TO_RAD, 0.0);
    double gamma = hinge.length();
    double phi;
    if (fabs(roll) > 0.000001)
        phi = atan(tilt/roll);
    else if (tilt > 0.0)
        phi = M_PI/2.0;
    else
        phi = -M_PI/2.0;

    if ((cos(phi) * roll < 0) || (sin(phi) * tilt < 0))
        gamma *= -1.0;


    // Compute Ti+1; equation 9
    matrix3x3 rz1 = rodrigues_formula(z, twist*DEG_TO_RAD/2.0 - phi);
    matrix3x3 ry2 = rodrigues_formula(y, gamma);
    matrix3x3 rz3 = rodrigues_formula(z, twist*DEG_TO_RAD/2.0 + phi);

    matrix3x3 tiplus1 = rz1*ry2;
    tiplus1 = tiplus1 * rz3;
    tiplus1 = tiplus1 * matrix3x3(x, y, z).transpose();

    // Compute Tmst; equation 10
    rz1 = rodrigues_formula(z, twist*DEG_TO_RAD/2.0 - phi);
    ry2 = rodrigues_formula(y, gamma/2.0);
    rz3 = rodrigues_formula(z, phi);

    matrix3x3 tmst = rz1*ry2;
    tmst = tmst * rz3;
    tmst = tmst * matrix3x3(x, y, z).transpose();

    // Compute the new origin; Equation 11
    vector3 new_origin = old_origin + shift * tmst.GetColumn(0) + slide * tmst.GetColumn(1) + rise * tmst.GetColumn(2);

    // return origin and direction vectors
    vector<vector3> results = {new_origin, tiplus1.GetColumn(0), tiplus1.GetColumn(1), tiplus1.GetColumn(2)};

    return results;
}


void HelicalParameters::ReferenceFrameToHelicalParameters(vector3 origin2, vector3 x2, vector3 y2, vector3 z2) {
    vector3 origin1(0.0, 0.0, 0.0);
    vector3 x1(1.0, 0.0, 0.0);
    vector3 y1(0.0, 1.0, 0.0);
    vector3 z1(0.0, 0.0, 1.0);

    matrix3x3 direction_vectors1 = matrix3x3(x1, y1, z1).transpose();
    matrix3x3 direction_vectors2 = matrix3x3(x2, y2, z2).transpose();

    // Compute the helical axis
    vector3 helical_axis = OpenBabel::cross(direction_vectors2.GetColumn(0) - direction_vectors1.GetColumn(0), direction_vectors2.GetColumn(1) - direction_vectors1.GetColumn(1));
    if (helical_axis.length() < 0.000001)
        helical_axis = vector3(0.0, 0.0, 1.0);
    else
        helical_axis /= helical_axis.length();

    // Compute TipInclination angle; Equation 17
    double Lambda1;
    if ((helical_axis.length() < 0.000001) || (direction_vectors1.GetColumn(2).length() < 0.000001))
        Lambda1 = 0.0;
    else {
        double dot = OpenBabel::dot(helical_axis/helical_axis.length(), direction_vectors1.GetColumn(2)/direction_vectors1.GetColumn(2).length());
        if (dot >= 1.0)
            Lambda1 = 0.0;
        else if (dot <= -1.0)
            Lambda1 = M_PI;
        else
            Lambda1 = acos(dot);
    }

    // Compute the TipInclination axis; Equation 18
    vector3 hinge1 = OpenBabel::cross(helical_axis, direction_vectors1.GetColumn(2));
    if (hinge1.length() > 0.000001)
        hinge1 /= hinge1.length();

    // Again compute TipInclination for the second base pair
    double Lambda2 = acos(OpenBabel::dot(helical_axis, direction_vectors2.GetColumn(2)));
    vector3 hinge2 = OpenBabel::cross(helical_axis, direction_vectors2.GetColumn(2));
    if (hinge2.length() > 0.000001)
        hinge2 /= hinge2.length();

    // Compute Ti'; Equation 19
    matrix3x3 direction_vectors1_h = rodrigues_formula(hinge1, -Lambda1) * direction_vectors1;
    matrix3x3 direction_vectors2_h = rodrigues_formula(hinge2, -Lambda2) * direction_vectors2;

    // Compute phi; Equation 20
    double phi_double_prime = acos(OpenBabel::dot(hinge1, direction_vectors1_h.GetColumn(1)));
    if (OpenBabel::dot(OpenBabel::cross(hinge1, direction_vectors1_h.GetColumn(1)), direction_vectors1_h.GetColumn(2)) > 0.0)
        phi_double_prime = fabs(phi_double_prime);
    else
        phi_double_prime = -fabs(phi_double_prime);

    // Compute tip and inclination; Equation 21
    tip = Lambda1 * cos(phi_double_prime) / DEG_TO_RAD;
    inclination = Lambda1 * sin(phi_double_prime) / DEG_TO_RAD;

    // Get axes for the twist
    vector3 t1 = direction_vectors1_h.GetColumn(1);
    t1 = t1 - OpenBabel::dot(t1, helical_axis) * helical_axis;
    t1 = t1 / t1.length();

    vector3 t2 = direction_vectors2_h.GetColumn(1);
    t2 = t2 - OpenBabel::dot(t2, helical_axis) * helical_axis;
    t2 = t2 / t1.length();

    // Compute twist; Equation 22
    h_twist = fabs(acos(OpenBabel::dot(t1, t2))) / DEG_TO_RAD;
    if (OpenBabel::dot(OpenBabel::cross(t1, t2), helical_axis) < 0.0)
        h_twist *= -1.0;

    // Compute translational parameters

    t2 = origin2 - origin1;
    // Compute rise; Equation 24
    h_rise = OpenBabel::dot(t2, helical_axis);

    t1 = t2 - h_rise * helical_axis;
    t1 = t1 - OpenBabel::dot(t1, helical_axis) * helical_axis;

    vector3 origin1_h;
    if (fabs(h_twist) < 0.000001)
        origin1_h = origin1 + 0.5 * t1;
    else {
        matrix3x3 temp = rodrigues_formula(helical_axis, M_PI/2.0 - h_twist*DEG_TO_RAD/2.0);
        vector3 axis = temp*t1;
        if (axis.length() > 0.0000001)
            axis /= axis.length();
        double s = 0.5 * t1.length() / sin(h_twist*DEG_TO_RAD / 2.0);
        origin1_h = origin1 + s * axis;
    }

    vector3 origin2_h = origin1_h + h_rise + helical_axis;

    // Compute x-displacement and y-displacement
    t1 = origin1_h - origin1;
    vector3 temp = OpenBabel::dot(t1, direction_vectors1_h.GetColumn(0));
    x_displacement = -1.0 * (temp.x() + temp.y() + temp.z());

    temp = OpenBabel::dot(t1, direction_vectors1_h.GetColumn(1));
    y_displacement = -1.0 * (temp.x() + temp.y() + temp.z());

    return;
}

void Backbone::deleteVectorAtom() {
    // Delete the extra atom that defines the vector connecting the backbone to the base
    // If it is already deleted, do nothing
    if (!vector_atom_deleted) {
        // Save the IDs of atoms connecting the two backbones and the ID of the atom connecting to the base
        array<size_t, 3> tmpIDs{backbone.GetAtom(interconnects[0])->GetId(),
                                          backbone.GetAtom(interconnects[1])->GetId(),
                                          backbone.GetAtom(linker[0])->GetId()};
        // Delete the atom
        backbone.DeleteAtom(getVector());
        // Get the new indices after the deletion of the atom
        interconnects = {backbone.GetAtomById(tmpIDs[0])->GetIdx(), backbone.GetAtomById(tmpIDs[1])->GetIdx()};
        linker = {backbone.GetAtomById(tmpIDs[2])->GetIdx(), 0};
        vector_atom_deleted = true;
    }
}

Backbone::Backbone(std::string file_path, std::array<unsigned, 2> interconnects, std::array<unsigned,2> linker, std::vector<std::vector<unsigned>> fixed_bonds) {
    // Constructor for the backbone

    vector_atom_deleted = false;
    this->interconnects = interconnects;
    this->linker = linker;
    this->fixed_bonds = fixed_bonds;

    // Read backbone file
    OBConversion conv(file_path);
    OBFormat *inFormat = conv.FormatFromExt(file_path);
    if (inFormat) {
        conv.SetInFormat(inFormat);
    } else {
        cerr << "Backbone: Cannot determine file type from extension. Defaulting to PDB." << endl;
        conv.SetInFormat("PDB");
    }
    if (!conv.Read(&backbone)) {
        throw std::runtime_error("Backbone: There was an error reading file for backbone: " + file_path);
        
    }

    // Basic validation of the indices
    validate();
}

void Backbone::validate() {

    // Make sure only two atoms are specified for connecting the backbones
    auto vec = interconnects;
    if (vec.size() != 2) {
        string error = "Backbone: Incorrect number of elements specified for the atoms connecting the two backbones.";
        error += " There must be exactly 2 indices indicated.";
        throw std::runtime_error(error);
    }
    copy(vec.begin(), vec.end(), interconnects.begin());

    // Make sure only two atoms are specified for connecting backbone to base
    vec = linker;
    if (vec.size() != 2) {
        string error = "Backbone: Incorrect number of elements specified for the atoms forming the vector connecting the backbone to the nucleobase.";
        error += " There must be exactly 2 indices indicated.";
        throw std::runtime_error(error);
    }
    copy(vec.begin(), vec.end(), linker.begin());

    for (auto i : interconnects)
        // Make sure the index is within the number of atoms in the backbone
        if (i > backbone.NumAtoms()) {
            throw std::runtime_error("Backbone: Interconnects atom indices are out of bounds.");
        }

    for (auto i : linker)
        // Make sure the index is within the number of atoms in the backbone
        if (i > backbone.NumAtoms()) {
            throw std::runtime_error("Backbone: Linker atom indices are out of bounds.");
        }

    OBAtom *link = backbone.GetAtom(static_cast<int>(linker[0])),
            *link2 = backbone.GetAtom(static_cast<int>(linker[1]));

    // Make sure the two provided atoms are neighbors
    if (!link->IsConnected(link2)) {
        throw std::runtime_error("Backbone: Atoms in defined in the Linker must be neighbors.");
    }

}


Base::Base(std::string nameT, std::string codeT, std::string file_path, std::array<std::size_t, 2> linkerT, std::string pairT) {
    // Constructor for the base

    name = nameT;
    code = codeT;
    linker = linkerT;

    if (file_path.empty()) {
        throw std::runtime_error("Bases: No file path was specified.");
    }

    // Read file specifying the geometry of the base
    OBConversion conv;
    OBFormat *in_format = conv.FormatFromExt(file_path);
    if (in_format)
        conv.SetInFormat(in_format);
    else {
        cerr << "Base: Cannot determine file type from extension. Defaulting to PDB." << endl;
        conv.SetInFormat("PDB");
    }
    OBMol molT;
    if (!conv.ReadFile(&molT, file_path)) {
        throw std::runtime_error("Base: There was an error reading file for base: " + file_path);
    }

    base = OBMol(molT);
    vector_atom_deleted = false;
    pair_name = pairT;

    // validate
    validate();
}


void Base::validate() {
    // Basic validation

    // Make sure we have the correct size for the atom indices
    if (linker.size() == 2) {
        for (auto i : linker)
            // Make sure the indices are within the number of atoms in the base
            if (i > base.NumAtoms()) {
                throw std::runtime_error("Base: Linker atom indices are out of bounds.");
            }

        OBAtom *link = base.GetAtom(static_cast<int>(linker[0])),
               *vec = base.GetAtom(static_cast<int>(linker[1]));

        // Make sure the atoms are neighbors
        if (!link->IsConnected(vec)) {
            throw std::runtime_error("Base: Atoms in defined in the Linker must be neighbors.");
        }

    } else {
        throw std::runtime_error("Base: Linker is of incorrect dimension.");
    }
}

Bases::Bases(std::vector<Base> input_bases) {

    // Loop through the bases and create a vector of them 
    for (unsigned i = 0; i < input_bases.size(); ++i) {
        string file_path = input_bases[i].file_path;

        auto vec = input_bases[i].linker;
        if (vec.size() != 2) {
            string error = "Incorrect number of elements specified for backbone connection. ";
            error += "There must be exactly 2 indices indicated.";
            throw std::runtime_error(error);
        }
        array<size_t, 2> linkers = {vec[0], vec[1]};

        auto name = input_bases[i].name;
        transform(name.begin(), name.end(), name.begin(), ::tolower);
        auto code = input_bases[i].code;
        transform(code.begin(), code.end(), code.begin(), ::toupper);
        auto pair_name = input_bases[i].pair_name;
        transform(pair_name.begin(), pair_name.end(), pair_name.begin(), ::tolower);

        // Call Base again to make sure we have the openbabel molecules
        bases.push_back(Base(name,code,file_path,linkers, pair_name));
    }

    // Checking to see if the bases have pairs
    all_bases_pair = true;
    for (auto v : bases) {
        auto pair_name = v.getBasePairName();
        transform(pair_name.begin(), pair_name.end(), pair_name.begin(), ::tolower);
        bool pair_found = false;
        for (auto w : bases) {
            if (w.getName().find(pair_name) != string::npos) {
                auto base_name = v.getName();
                transform(base_name.begin(), base_name.end(), base_name.begin(), ::tolower);
                name_base_map.insert(pair<string,Base>(v.getName(),w));
                pair_found = true;
                break;
            }
        }
        if (!pair_found) {
            all_bases_pair = false;
            break;
        }
    }
}

std::vector<Base> Bases::getBasesFromStrand(std::vector<std::string> strand) {
    // Return the bases given the names of the bases in the strand
    vector<Base> bases;
    for (auto &v : strand) {
        transform(v.begin(), v.end(), v.begin(), ::tolower);
        bases.push_back(getBaseFromName(v));
    }
    return bases;
}

std::vector<Base> Bases::getComplimentBasesFromStrand(std::vector<std::string> strand) {
    // Return the complimentary bases given the names of the bases in the strand
    if (all_bases_pair) {
        vector<Base> compliment_bases;
        for (auto &v : strand) {
            transform(v.begin(), v.end(), v.begin(), ::tolower);
            compliment_bases.push_back(name_base_map.find(v)->second);
        }
        return compliment_bases;
    }

    // If not all bases have complimentary bases, return empty vector
    return vector<Base>();
}

BaseUnit::BaseUnit(Base base, Backbone backbone, double glycosidic_bond_distance) {

    // We go ahead and center everything so we can rotate later
    backbone.center();

    // Setting up an array of necessary OBAtom for alignment
    array< OBAtom*, 4 > atoms{base.getLinker() , base.getVector()  ,
                                   backbone.getLinker(), backbone.getVector() };

    vector<vector3> ref {atoms[1]->GetVector(), atoms[0]->GetVector()},
            target {atoms[2]->GetVector(), atoms[3]->GetVector()},
            result;

    // Aligning backbone linker to base linker
    OBAlign align(ref, target);
    matrix3x3 matrix;
    if (align.Align()) {
        matrix = align.GetRotMatrix();
    } else {
        throw std::runtime_error("Failed to align backbone to base. Check base and backbone.");
    }
    double rot[9];
    matrix.GetArray(rot);

    backbone.rotate(rot);

    // Necessary to invert to get correct DNA and RNA structures
    // if alignment requires improper rotation
    if (matrix.determinant() < 0) {
        double invert[9] = {1, 0, 0, 0, 1, 0, 0, 0, -1};
        backbone.rotate(invert);
    }

    // The bond between the base and the backbone is determined by the
    // distance between the two base linkers
    backbone.translate(atoms[1]->GetVector() - atoms[2]->GetVector());

    // Get indices of fixed atoms
    unsigned num_atoms = base.getMolecule().NumAtoms() - 1;
    for (int i=0; i < backbone.fixed_bonds.size(); i++) {
        fixed_bonds.push_back({});
        for (auto j: backbone.fixed_bonds[i]) {
            if (j == backbone.linker[1]) {
                // This atom will be deleted, so put the index of the atom from the base
                if (base.linker[0] < base.linker[1])
                    fixed_bonds[i].push_back(base.linker[0]);
                else
                    fixed_bonds[i].push_back(base.linker[0]-1);
            }
            else if (j < backbone.linker[1])
                // The index of the atom will not change because of deletion
                fixed_bonds[i].push_back(j + num_atoms);
            else
                // Subtract one because of the deletion of the atom
                fixed_bonds[i].push_back(j-1 + num_atoms);
        }
    }

    // Deleting old atoms that formed the vector
    backbone.deleteVectorAtom();
    base.deleteVectorAtom();

    // Create the molecule
    OBMol mol;
    mol += base.getMolecule();

    num_atoms = mol.NumAtoms();
    mol += backbone.getMolecule();

    // Form a bond between base and backbone
    mol.AddBond(base.getLinker()->GetIdx(), backbone.getLinker()->GetIdx() + num_atoms, 1);
    OBBond* new_bond = mol.GetBond(base.getLinker()->GetIdx(), backbone.getLinker()->GetIdx() + num_atoms);
    // Set Length prints useless messages to the console. We need to supress it
    OBMessageHandler handler = OBMessageHandler();
    handler.StartErrorWrap();
    if (glycosidic_bond_distance != 0.0)
        new_bond->SetLength(glycosidic_bond_distance);
    else
        new_bond->SetLength(mol.GetAtom(base.getLinker()->GetIdx()), new_bond->GetEquibLength());
    handler.StopErrorWrap();

    // garbage code for debugging
    // TODO delete this code
    //{
    //    OBConversion conv;
    //    std::filebuf fb;
    //    fb.open("base_unit_example_out.pdb", std::ios::out);
    //    std::ostream fileStream(&fb);
    //    conv.SetOutFormat("PDB");
    //    conv.SetOutStream(&fileStream);
    //    OBMol mol;
    //    mol += base.getMolecule();
    //    mol += backbone.getMolecule();
    //    mol.AddBond(base.getLinker()->GetIdx(), backbone.getLinker()->GetIdx() + num_atoms, 1);
    //    OBBond* new_bond = mol.GetBond(base.getLinker()->GetIdx(), backbone.getLinker()->GetIdx() + num_atoms);
    //    new_bond->SetLength(mol.GetAtom(base.getLinker()->GetIdx()), new_bond->GetEquibLength());
    //    conv.Write(&mol);
    //}

    // Set variables
    unit = mol;
    base_connect_index = num_atoms + backbone.getLinker()->GetIdx();
    backbone_interconnects = {backbone.getHead()->GetIdx() + num_atoms, backbone.getTail()->GetIdx() + num_atoms};
    base_index_range = {1, num_atoms};
    backbone_index_range = {num_atoms + 1, mol.NumAtoms()};
}
