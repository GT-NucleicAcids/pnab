//
// Created by jbarnett8 on 9/20/17.
//

#ifndef PNAB_FILEPARSER_H
#define PNAB_FILEPARSER_H

#include <string>
#include <vector>
#include <array>
#include <iostream>
#include <sstream>
#include <map>
#include <algorithm>

/**
 * \brief Used to hold the data for a generic field within a Category
 * @tparam T Any typename
 */
template <typename T>
class Field {
private:
    std::string name;                                           //!< The name used to identify the Field
    T data{};                                                     //!< A generic container for the data of the Field
    bool value_set = false;
    bool is_required = true;

public:

    /**
     * \brief Creates a Field object with a given name
     * @param name The intended name for the Field
     */
    Field(std::string name, bool is_required = true) {
        this->is_required = is_required;
        setName(name);
    };

    /**
     * \brief Changes name of the Field
     * @param str New name for the Field
     */
    void setName(std::string str) {
        std::transform(name.begin(), name.end(), name.begin(), ::tolower);
        name = str;
    };

    /**
     * \brief Gives the name of the field
     * @return The name of the Field
     */
    std::string getName() {
        return name;
    };

    /**
     * \brief Sets the data for the Field
     * @param d The data to give the Field
     */
    void setData(T d) {
        data = d;
        value_set = true;
    };

    /**
     * \brief Gives the data the Field holds
     * @return The data of the Field
     */
    T getData() {
        return data;
    };

    /**
     * \brief Returns whether or not the field has been set
     * @return True if a value has been set, false otherwise
     */
    bool isSet() {
        return value_set;
    }

    /**
     * \brief Set where the field is required or not.
     * @param val true if required, false if not
     */
    void setIsRequired(bool val) {
        is_required = val;
    }

    /**
     * \brief Returns whether or not the field is required.
     * @return true is the Field is required, false otherwise
     */
    bool isRequired() {
        return is_required;
    }

    /**
     * \brief Takes in a string with all values to the right of the "=" in the input file which is then interpreted
     * and stored in data with a function that depends on the data type.
     * @param str The string containing all values on the right-hand-size of the "=" from an input file
     */
    virtual void parse(std::string str) = 0;
};

/**
 * \brief Used to hold a positive integer value
 */
class SizeField : public Field< std::size_t > {
public:
    using Field::Field;
    void parse(std::string str) {
        while (!isalnum(str.front()))
            str.erase(0, 1);
        if (str.empty()) {
            std::cerr << "Field must contain alpha-numeric characters" << std::endl;
            exit(1);
        }
        setData(static_cast<std::size_t>( std::stoul( str ) ) );
    };
};

/**
 * \brief Used to hold multiple positive integer values
 */
class SizeVecField : public Field< std::vector< std::size_t > > {
public:
    using Field::Field;
    void parse(std::string str) {
        std::vector< std::size_t > vec;
        std::istringstream ss(str);
        std::string token;
        while ( getline(ss, token, ',') ) {
            while (!isalnum(token.front()) && token.front() != '.')
                token.erase(0, 1);
            if (token.empty()) {
                std::cerr << "Field must contain alpha-numeric characters. " << std::endl;
                exit(1);
            }
            vec.push_back( std::stoul(token) );
        }
        setData( vec );
    };
};

/**
 * \brief Used to hold a string value
 */
class StringField : public Field< std::string > {
public:
    using Field::Field;
    void parse(std::string str) {
        while (!isalnum(str.front()) && !isgraph(str.front()))
            str.erase(0, 1);
        if (str.empty()) {
            std::cerr << "Field must contain alpha-numeric characters" << std::endl;
            exit(1);
        }
        setData( str );
    };
};

/**
 * \brief Used to hold multiple string values
 */
class StringVecField : public Field< std::vector< std::string > > {
public:
    using Field::Field;
    void parse(std::string str) {
        std::vector< std::string > vec;
        std::istringstream ss(str);
        std::string token;
        while ( getline(ss, token, ',') ) {
            while (!isalnum(token.front()) && !isgraph(str.front()))
                token.erase(0, 1);
            if (token.empty()) {
                std::cerr << "Field must contain alpha-numeric characters. " << std::endl;
                exit(1);
            }
            vec.push_back( token );
        }
        setData( vec );
    };
};

/**
 * \brief Used to hold real number values
 */
class DoubleField : public Field< double > {
public:
    using Field::Field;
    void parse(std::string str) {
        auto front = str.front();
        while (!str.empty() && !isalnum(front) && !(front == '.' || front == '-'))
            str.erase(0, 1);
        if (str.empty()) {
            std::cerr << "Field must contain alpha-numeric characters. " << std::endl;
            exit(1);
        }
        setData( std::stod(str) );
    };
};

/**
 * \brief Used to hold multiple real number values
 */
class DoubleVecField : public Field< std::vector<double> > {
public:
    using Field::Field;
    void parse(std::string str) {
        std::vector< double > vec;
        std::istringstream ss(str);
        std::string token;
        while ( getline(ss, token, ',') ) {
            auto front = token.front();
            while (!token.empty() && !isalnum(front) && !(front == '.' || front == '-')) {
                token.erase(0, 1);
                front = token.front();
            }
            if (token.empty()) {
                std::cerr << "Field must contain alpha-numeric characters. " << std::endl;
                exit(1);
            }
            vec.push_back( std::stod(token) );
        }
        setData( vec );
    };
};

/**
 * \brief Used to define a single category and all fields within
 *
 * This class is designed to allow for the easy creation of new categories and easy registration of fields
 * within. Each type of field has its own Field type derived from the virtual Field class. If additional
 * Field types are required, they must be manually added as a new class derived from Field with a corresponding
 * addition to the FieldType enum.
 */
class Category {

public:

    /**
     * \brief Creates a \code{Category} with a given name
     * @param str The intended name for the category used as a section header for an input file
     */
    Category(std::string str) {
        std::transform(str.begin(), str.end(), str.begin(), ::toupper);
        name = str;
    };

    /**
     * Returns of the name of the category
     * @return
     */
    std::string getName() {
        return name;
    };

    /**
     * \brief Registers field into vector holding fields for positive integer values
     * @param name
     */
    void registerSizeField(std::string name, bool is_required = true) {
        std::transform(name.begin(), name.end(), name.begin(), ::tolower);
        size_fields.push_back(SizeField(name, is_required));
        stringToFieldTypeMap.insert(std::pair<std::string, FieldType>(name, SIZE_FIELD));
    };

    /**
     * \brief Registers field into vector holding fields for multiple positive integer values
     * @param name
     */
    void registerSizeVecField(std::string name, bool is_required = true) {
        std::transform(name.begin(), name.end(), name.begin(), ::tolower);
        size_vec_fields.push_back(SizeVecField(name, is_required));
        stringToFieldTypeMap.insert(std::pair<std::string, FieldType>(name, SIZE_VEC_FIELD));
    };

    /**
     * \brief Registers field into vector holding fields for string values
     * @param name
     */
    void registerStringField(std::string name, bool is_required = true) {
        std::transform(name.begin(), name.end(), name.begin(), ::tolower);
        string_fields.push_back(StringField(name, is_required));
        stringToFieldTypeMap.insert(std::pair<std::string, FieldType>(name, STRING_FIELD));
    };

    /**
     * \brief Registers field into vector holding fields for multiple string values
     * @param name
     */
    void registerStringVecField(std::string name, bool is_required = true) {
        std::transform(name.begin(), name.end(), name.begin(), ::tolower);
        string_vec_fields.push_back(StringVecField(name, is_required));
        stringToFieldTypeMap.insert(std::pair<std::string, FieldType>(name, STRING_VEC_FIELD));
    };

    /**
     * \brief Registers field into vector holding fields for real number values
     * @param name
     */
    void registerDoubleField(std::string name, bool is_required = true) {
        std::transform(name.begin(), name.end(), name.begin(), ::tolower);
        double_fields.push_back(DoubleField(name, is_required));
        stringToFieldTypeMap.insert(std::pair<std::string, FieldType>(name, DOUBLE_FIELD));
    };

    /**
     * \brief Registers field into vector holding fields for multiple real number values
     * @param name
     */
    void registerDoubleVecField(std::string name, bool is_required = true) {
        std::transform(name.begin(), name.end(), name.begin(), ::tolower);
        double_vec_fields.push_back(DoubleVecField(name, is_required));
        stringToFieldTypeMap.insert(std::pair<std::string, FieldType>(name, DOUBLE_VEC_FIELD));
    };

    /**
     * \brief Processes a full line of an input file after it has been stripped of all comments
     * @param line The line of the input file to process
     */
    void parseLine(std::string line);

    /**
     * \brief Checks whether all required fields are set within a Category.
     */
    void validate();

    /**
     * \brief Let us know if the category is empty or not
     * @return indicates if Category is empty or not
     */
    bool empty() {
        name.length() ? 0: false, true;
    }

    /**
     * \brief Prints all of the information pertaining to this Category
     */
    void printCategory();

    /**
     * \brief Returns the data within the field specified in the name
     * @param name The name of the field you want in the Category
     * @param default_val If the value in the field is not set, set to this value
     * @return The value, either actaul value or default
     */
    std::size_t getSizeFieldDataFromName(std::string name, std::size_t default_val = 0) {
        checkIfNameExists(name);
        auto it = findNameInFields(name, size_fields);
        if (it != size_fields.end())
            return it->isSet() ? it->getData() : default_val;
    }

    std::vector<std::size_t> getSizeVecFieldDataFromName(std::string name, std::vector<std::size_t> default_val = {}) {
        checkIfNameExists(name);
        auto it = findNameInFields(name, size_vec_fields);
        if (it != size_vec_fields.end())
            return it->isSet() ? it->getData() : default_val;
    }

    std::string getStringFieldDataFromName(std::string name, std::string default_val = "") {
        checkIfNameExists(name);
        auto it = findNameInFields(name, string_fields);
        if (it != string_fields.end())
            return it->isSet() ? it->getData() : default_val;
    }

    std::vector<std::string> getStringVecFieldDataFromName(std::string name, std::vector<std::string> default_val = {}) {
        checkIfNameExists(name);
        auto it = findNameInFields(name, string_vec_fields);
        if (it != string_vec_fields.end())
            return it->isSet() ? it->getData() : default_val;
    }

    double getDoubleFieldDataFromName(std::string name, double default_val = 0) {
        checkIfNameExists(name);
        auto it = findNameInFields(name, double_fields);
        if (it != double_fields.end())
            return it->isSet() ? it->getData() : default_val;
    }

    std::vector<double> getDoubleVecFieldDataFromName(std::string name, std::vector<double> default_val = {}) {
        checkIfNameExists(name);
        auto it = findNameInFields(name, double_vec_fields);
        if (it != double_vec_fields.end())
            return it->isSet() ? it->getData() : default_val;
    }

private:
    std::string name;
    std::vector< SizeField > size_fields;                           //!< Used to hold single positive integer value
    std::vector< SizeVecField > size_vec_fields;                    //!< Used to hold multiple positive integer values
    std::vector< StringField > string_fields;                       //!< Used to hold single string value
    std::vector< StringVecField > string_vec_fields;                //!< Used to hold multiple string values (i.e. chain of nucleoties)
    std::vector< DoubleField > double_fields;                       //!< Used to hold single real number value
    std::vector< DoubleVecField > double_vec_fields;                //!< Used to hold multiple real number values

    /**
     * \brief Returns an iterator pointed at the position where the name occurred within the vector of fields
     * @tparam T Any derived \code{Field} class
     * @param name The name of the field we are searching for within the vector of fields
     * @param fields The vector of fields
     * @return An iterator pointed at the position where the name occurred within \code{fields}
     */
    template <typename T>
    typename std::vector<T>::iterator findNameInFields(std::string name, std::vector<T> &fields) {
        auto p =  std::find_if(fields.begin(), fields.end(),
                            [name](T & t) -> bool {
                                return t.getName().compare(name) == 0;
                            }
        );
        if (p != fields.end())
            return p;
        else {
            std::cerr << "findNameInFields(): Name '" << name << "' does not exist in the field" << std::endl;
        }
    };

    void checkIfNameExists(std::string &name) {
        std::transform(name.begin(), name.end(), name.begin(), ::tolower);
        auto field_type_it = stringToFieldTypeMap.find(name);
        if (field_type_it == stringToFieldTypeMap.end()) {
            std::cerr << "Field \"" << name << "\" in category \"" << getName() << " is not registered."
                 << " Please check input file." << std::endl;
            exit(1);
        }
    }

    /**
     * \brief Used to distinguish between different field types in \code{stringToFieldTypeMap}
     */
    enum FieldType {
        SIZE_FIELD,
        SIZE_VEC_FIELD,
        STRING_FIELD,
        STRING_VEC_FIELD,
        DOUBLE_FIELD,
        DOUBLE_VEC_FIELD
    };
    std::map< std::string, FieldType > stringToFieldTypeMap;        //!< Used to find what field type a field within a category is

};

class FileParser {

public:
    /**
     * \brief Empty constructor. Populates available Category objects and Field objects within with available options
     */
    FileParser();

    /**
     * Same as the empty constructor, but provides file path to input file
     * @param path Path to input file
     */
    FileParser(std::string path) : FileParser() {
        file_path = path;
    };

    /**
     * \brief Reads the file given in the path \code{file_path} and stores data
     */
    void readFile();

    /**
     * \brief Sets the file path for the FileParser object used when calling readFile()
     * @param path
     */
    void setFilePath(std::string path) {
        file_path = path;
    };

    /**
     * \brief Print all available information contained within FileParser
     */
    void print() {
        for (auto cat : stringToCategoryMap)
            cat.second.printCategory();
    }

    /**
     * \brief Returns pointer to the Category like the private function getCategoryByName, but it's constant
     * @param name The name of the category you want
     * @return The category with name \code{name}
     */
    const Category& getCategory(std::string name) {
        return getCategoryByName(name);
    }

    /**
     * \brief Gives the number of base categories declared
     * @return The count of base categories
     */
    unsigned getNumBaseCategories() {
        return num_base_categories;
    }

private:
    /**
     * \brief Add a category with a given name which defines Category header in the input file
     * @param str
     */
    void registerCategory(Category cat) {
        stringToCategoryMap.insert(std::pair<std::string, Category>(cat.getName(), cat));
    }

    /**
     * \brief Get a category given a name
     * @param name Name of the category to search for
     * @return Pointer to the category we want
     */
    Category& getCategoryByName(std::string name) {
        std::transform(name.begin(), name.end(), name.begin(), ::toupper);
        auto it = stringToCategoryMap.find(name);
        if (it != stringToCategoryMap.end())
            return it->second;
        std::cerr << "getCategoryByName(): Category \"" << name << "\" does not exist." << std::endl;
        return stringToCategoryMap.end()->second;
    }

    std::map< std::string, Category > stringToCategoryMap;        //!< Used to hold list of Category class objects
    std::string file_path;                                        //!< String containing path to input file
    unsigned num_base_categories = 0;                             //!< Keeps track of the number of base parameters used

};

#endif //PNAB_FILEPARSER_H
