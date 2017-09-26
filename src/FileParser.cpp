//
// Created by jbarnett8 on 9/20/17.
//

#include <algorithm>
#include "FileParser.h"

using namespace std;

void FileParser::FileParser() {
    Category runtimeParameters("RUNTIME PARAMETERS");
    // Geometric properties of bases
    vector<string> single_doubles_geometric{"Rise", "X_Disp, Y_Disp, Inclination, Tip, Twist"};
    for (auto str : single_doubles_geometric)
        runtimeParameters.registerDoubleField(str);

    // Energetic filter properties of conformations
    vector<string> single_doubles_energy{"Max_Total_Energy", "Max_Angle_Energy", "Max_Bond_Energy",
                                         "Max_VDW_Energy", "Max_Torsion_Energy", "Max_Distance"};
    for (auto str : single_doubles_energy)
        runtimeParameters.registerDoubleField(str);

    // Force-Field parameters
    vector<string> single_strings_ff{"Type", "Parameter_File"};
    for (auto str : single_strings_ff)
        runtimeParameters.registerStringField(str);
    runtimeParameters.registerDoubleField("Base_to_Backbone_Bond_Length")

    // Search Algorithm parameters
    runtimeParameters.registerStringField("Algorithm");
    vector<string> single_size_search{"Search_Size", "Step_Size", };

    registerCategory(runtimeParameters);
}

void Category::parseLine(std::string line) {
    string field = line.substr(0, line.find("="));
    string value = line.substr(line.find("=") + 1);

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
            auto it = findNameInFields(name, size_fields);
            if (it != size_fields.end())
                it->parse(value);
            break;
        }
        case SIZE_VEC_FIELD: {
            auto it = findNameInFields(name, size_vec_fields);
            if (it != size_vec_fields.end())
                it->parse(value);
            break;
        }
        case STRING_FIELD: {
            auto it = findNameInFields(name, string_fields);
            if (it != string_fields.end())
                it->parse(value);
            break;
        }
        case STRING_VEC_FIELD: {
            auto it = findNameInFields(name, string_vec_fields);
            if (it != string_vec_fields.end())
                it->parse(value);
            break;
        }
        case DOUBLE_FIELD: {
            auto it = findNameInFields(name, double_fields);
            if (it != double_fields.end())
                it->parse(value);
            break;
        }
        case DOUBLE_VEC_FIELD: {
            auto it = findNameInFields(name, double_vec_fields);
            if (it != double_vec_fields.end())
                it->parse(value);
            break;
        }
    }
}

void Category::printCategory() {
    std::size_t printLength = name.length() + std::string("|  Printing Category: ").length();
    cout << std::string("-", printLength + 3) << endl;
    cout << "|  Printing Category: " << name << "  |" << endl;
    cout << std::string("-", printLength + 3) << endl;
    for (auto f : stringToFieldTypeMap) {
        switch (f.second) {
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
