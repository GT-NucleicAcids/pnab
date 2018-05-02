//
// Created by jbarnett8 on 9/20/17.
//

#include <algorithm>
#include <fstream>
#include "FileParser.h"

using namespace std;

FileParser::FileParser() {

    // Geometric properties of bases //////////////////////////////////////////////////////
    Category runtimeParameters("RUNTIME PARAMETERS");

    // Energetic filter properties of conformations
    vector<string> single_doubles_energy{"Max_Total_Energy", "Max_Angle_Energy", "Max_Bond_Energy",
                                         "Max_VDW_Energy", "Max_Torsion_Energy", "Max_Backbone_Interlink_Distance"};
    for (auto str : single_doubles_energy)
        runtimeParameters.registerDoubleField(str);

    // Force-Field parameters
    vector<string> single_strings_ff{"Force_Field_Type", "Force_Field_Parameter_File"};
    for (auto str : single_strings_ff)
        runtimeParameters.registerStringField(str, false);
    runtimeParameters.registerDoubleField("Base_to_Backbone_Bond_Length", false);

    // Search Algorithm parameters
    runtimeParameters.registerStringField("Algorithm");
    vector<string> single_size_search{"Search_Size", "Dihedral_Step_Size", "Search_Step_Size",
                                      "Chain_Length"};
    vector<bool> is_needed{true, false, false, false};
    size_t i = 0;
    for (auto str : single_size_search) {
        runtimeParameters.registerSizeField(str, is_needed[i++]);
    }
    runtimeParameters.registerStringVecField("Strand_Base_Names", false);
    runtimeParameters.registerStringField("Is_Double_Stranded", false);

    // Store category in FileParser
    registerCategory(runtimeParameters);

    // Helical Parameters //////////////////////////////////////////////////////////////////
    Category helical_parameters("HELICAL PARAMETERS");
    vector<string> vec_doubles_hp{"Tilt", "Roll", "Twist", "Shift", "Slide", "Rise", "Buckle", "Propeller", "Opening",
                                  "Shear", "Stretch", "Stagger", "Inclination", "Tip", "X_Displacement",
                                  "Y_Displacement"};
    for (auto v : vec_doubles_hp)
        helical_parameters.registerDoubleVecField(v, v.compare("Twist") == 0 || v.compare("Rise") == 0);
    registerCategory(helical_parameters);

    // Backbone Parameters /////////////////////////////////////////////////////////////////
    Category backboneParameters("BACKBONE PARAMETERS");
    vector<string> mult_size_backbone{"Interconnects", "Base_Connect"};
    for (auto str : mult_size_backbone)
        backboneParameters.registerSizeVecField(str);
    backboneParameters.registerStringField("Backbone_File_Path");

    registerCategory(backboneParameters);

/**********************************************************************************************
//    Base parameters are added as needed in readFile()
//    // Base Parameters /////////////////////////////////////////////////////////////////////
//    Category base0Parameters("BASE PARAMETERS");
//    vector<string> mult_string_base{"Code", "Name", "Base_File_Path"};
//    for (auto str : mult_string_base)
//        base0Parameters.registerStringField(str);
//    base0Parameters.registerSizeVecField("Backbone_Connect");

//    registerCategory(base0Parameters);
***********************************************************************************************/
}

void FileParser::readFile() {
    map<string,Category> found_categories;
    ifstream f(file_path);
    if (f.is_open()) {
        cout << "Reading: " << file_path << endl;
        string line;
        size_t line_number = 0;
        Category *c;
        while (getline(f, line)) {
            line_number++;
            size_t comment_start = line.find_first_of("#");

            // Can't have length of -1 in substr(pos,len)
            if (comment_start == 0)
                continue;
            string no_comment = line.substr(0, line.find_first_of("#"));

            // If we can't find anything besides whitespace, we don't care about that line
            if (no_comment.find_first_not_of(" \t\n") == string::npos)
                continue;

            // Strip trailing whitespace
            no_comment = no_comment.substr(0, no_comment.find_last_not_of(" \t\n") + 1);
            // Strip leading whitespace
            no_comment = no_comment.substr(no_comment.find_first_not_of(" \t\n"));

            if (no_comment.find("=") != string::npos && !c->empty()) {
                c->parseLine(no_comment);
            } else if (no_comment.find("=") != string::npos && c->empty()) {
                cerr << "Declared field before declaring a category on line " << line_number
                     << "with text \"" << line << "\". Please place field under appropriate category"
                     << endl;
                exit(1);
            }
            else {
                transform(no_comment.begin(),no_comment.end(),no_comment.begin(),::toupper);
                auto it = stringToCategoryMap.find(no_comment);
                if (found_categories.find(no_comment) != found_categories.end()
                    || no_comment.find("BASE PARAMETERS") != string::npos) {
                    // I find this workaround unsavory. The input file is easier to use if we can
                    // simply have multiple instances of the BASE PARAMETERS category, so I'm hacking my way
                    // through to make it possible. An more ideal solution would involve refactoring the FileParser
                    // class entirely.
                    if (no_comment.find("BASE PARAMETERS") == string::npos) {
                        cerr << "Category \"" << no_comment << "\" on line " << line_number
                             << " was already used. Please check the input file." << endl;
                        exit(1);
                    }

                    // Base Parameters /////////////////////////////////////////////////////////////////////
                    num_base_categories++;
                    Category baseParameters("BASE PARAMETERS " + to_string(num_base_categories));
                    vector<string> mult_string_base{"Code", "Name", "Base_File_Path", "Pair_Name"};
                    for (auto str : mult_string_base)
                        baseParameters.registerStringField(str);
                    baseParameters.registerSizeVecField("Backbone_Connect");

                    registerCategory(baseParameters);
                    it = stringToCategoryMap.find("BASE PARAMETERS " + to_string(num_base_categories));
                }
                if (it == stringToCategoryMap.end()) {
                    cerr << "Category \"" << no_comment << "\" on line " << line_number
                         << " does not exist. Please check the input file." << endl;
                    exit(1);
                }
                c = &it->second;
                found_categories.insert(*it);
            }
        }
        for (auto category_map : stringToCategoryMap) {
            category_map.second.validate();
        }
    } else {
        cerr << "There was an error opening file: " << file_path << endl;
        exit(1);
    }
}

void Category::parseLine(std::string line) {
    string field = line.substr(0, line.find("="));
    string value = line.substr(line.find("=") + 1);

    transform(field.begin(), field.end(), field.begin(), ::tolower);

    if (value.find_first_not_of(" \n\t") == string::npos) {
        cerr << "Error: Empty field \"" << field << "\" in category \"" << getName() << "." << endl;
        exit(1);
    }

    // Need to make sure the field actually exists within the map somewhere. A nice benefit to this is that if it
    // exists in the map, it is guaranteed to exist with the vector of Fields as well unless the functions of
    // std::vector fail, in which case there are bigger problems...
    auto field_type_it = stringToFieldTypeMap.find(field);
    if (field_type_it == stringToFieldTypeMap.end()) {
        cerr << "Field \"" << field << "\" in category \"" << getName() << " is not registered."
             << " Please check input file." << endl;
        exit(1);
    }
    FieldType type = field_type_it->second;
    switch (type) {
        case SIZE_FIELD: {
            auto it = findNameInFields(field, size_fields);
            if (it != size_fields.end())
                it->parse(value);
            break;
        }
        case SIZE_VEC_FIELD: {
            auto it = findNameInFields(field, size_vec_fields);
            if (it != size_vec_fields.end())
                it->parse(value);
            break;
        }
        case STRING_FIELD: {
            auto it = findNameInFields(field, string_fields);
            if (it != string_fields.end())
                it->parse(value);
            break;
        }
        case STRING_VEC_FIELD: {
            auto it = findNameInFields(field, string_vec_fields);
            if (it != string_vec_fields.end())
                it->parse(value);
            break;
        }
        case DOUBLE_FIELD: {
            auto it = findNameInFields(field, double_fields);
            if (it != double_fields.end())
                it->parse(value);
            break;
        }
        case DOUBLE_VEC_FIELD: {
            auto it = findNameInFields(field, double_vec_fields);
            if (it != double_vec_fields.end())
                it->parse(value);
            break;
        }
    }
}

void Category::printCategory() {
    std::size_t printLength = name.length() + std::string("|  Printing Category: ").length();
    cout << std::string(printLength + 3, '-') << endl;
    cout << "|  Printing Category: " << name << "  |" << endl;
    cout << std::string(printLength + 3, '-') << endl;
    std::vector<FieldType> types{SIZE_FIELD, SIZE_VEC_FIELD, STRING_FIELD,
                                 STRING_VEC_FIELD, DOUBLE_FIELD,
                                 DOUBLE_VEC_FIELD};
    for (auto f : types) {
        switch (f) {
            case SIZE_FIELD: {
                for (auto sf : size_fields) {
                    cout << "\t" << sf.getName() << ": " << sf.getData() << endl;
                }
                break;
            }
            case SIZE_VEC_FIELD: {
                for (auto sf : size_vec_fields) {
                    cout << "\t" << sf.getName() << ": ";
                    for (auto itemS : sf.getData())
                        cout << itemS << ", ";
                    cout << endl;
                }
                break;
            }
            case STRING_FIELD: {
                for (auto sf : string_fields) {
                    cout << "\t" << sf.getName() << ": " << sf.getData() << endl;
                }
                break;
            }
            case STRING_VEC_FIELD: {
                for (auto sf : string_vec_fields) {
                    cout << "\t" << sf.getName() << ": ";
                    for (auto itemStr : sf.getData())
                        cout << itemStr << ", ";
                    cout << endl;
                }
                break;
            }
            case DOUBLE_FIELD: {
                for (auto sf : double_fields) {
                    cout << "\t" << sf.getName() << ": " << sf.getData() << endl;
                }
                break;
            }
            case DOUBLE_VEC_FIELD: {
                for (auto sf : double_vec_fields) {
                    cout << "\t" << sf.getName() << ": ";
                    for (auto itemD : sf.getData())
                        cout << itemD << ", ";
                    cout << endl;
                }
                break;
            }
        }
    }
    cout << endl;
}

void Category::validate() {
    std::vector<FieldType> types{SIZE_FIELD, SIZE_VEC_FIELD, STRING_FIELD,
                                 STRING_VEC_FIELD, DOUBLE_FIELD,
                                 DOUBLE_VEC_FIELD};
    bool error_found = false;
    for (auto f : types) {
        switch (f) {
            case SIZE_FIELD: {
                for (auto sf : size_fields) {
                    if (!sf.isSet() & sf.isRequired()) {
                        cout << "Required field \"" << sf.getName() << "\" in Category \"" << name
                             << "\" is not set." << endl;
                        error_found = true;
                    }
                }
                break;
            }
            case SIZE_VEC_FIELD: {
                for (auto sf : size_vec_fields) {
                    if (!sf.isSet() & sf.isRequired()) {
                        cout << "Required field \"" << sf.getName() << "\" in Category \"" << name
                             << "\" is not set." << endl;
                        error_found = true;
                    }
                }
                break;
            }
            case STRING_FIELD: {
                for (auto sf : string_fields) {
                    if (!sf.isSet() & sf.isRequired()) {
                        cout << "Required field \"" << sf.getName() << "\" in Category \"" << name
                             << "\" is not set." << endl;
                        error_found = true;
                    }
                }
                break;
            }
            case STRING_VEC_FIELD: {
                for (auto sf : string_vec_fields) {
                    if (!sf.isSet() & sf.isRequired()) {
                        cout << "Required field \"" << sf.getName() << "\" in Category \"" << name
                             << "\" is not set." << endl;
                        error_found = true;
                    }
                }
                break;
            }
            case DOUBLE_FIELD: {
                for (auto sf : double_fields) {
                    if (!sf.isSet() & sf.isRequired()) {
                        cout << "Required field \"" << sf.getName() << "\" in Category \"" << name
                             << "\" is not set." << endl;
                        error_found = true;
                    }
                }
                break;
            }
            case DOUBLE_VEC_FIELD: {
                for (auto sf : double_vec_fields) {
                    if (!sf.isSet() & sf.isRequired()) {
                        cout << "Required field \"" << sf.getName() << "\" in Category \"" << name
                             << "\" is not set." << endl;
                        error_found = true;
                    }
                }
                break;
            }
        }
    }
    // We want to be able to let user know all of the fields that are missing all at once, so we exit here
    if (error_found)
        exit(1);
}