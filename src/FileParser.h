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
    FileParser(std::string path);

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

private:
    /**
     * \brief Add a category with a given name which defines Category header in the input file
     * @param str
     */
    void registerCategory(Category cat) {
        std::transform(name.begin(), name.end(), name.begin(), ::toupper);
        stringToCategoryMap.insert(std::pair<std::string, Category>(cat.getName(), cat));
    }

    /**
     * \brief Gives an iterator pointing to the matching element in the map
     * @param name The name of the category to search for
     * @return An iterator pointing to the matching element in the map
     */
    std::map< std::string, Category >::iterator getCategoryFromName(std::string name) {
        return stringToCategoryMap.find(name);
    }


private:
    std::map< std::string, Category > stringToCategoryMap;        //!< Used to hold list of Category class objects
    std::string file_path;                                        //!< String containing path to input file

};

/**
 * \brief Used to hold the data for a generic field within a Category
 * @tparam T Any typename
 */
template <typename T>
class Field {
private:
    std::string name;                                           //!< The name used to identify the Field
    T data;                                                     //!< A generic container for the data of the Field
public:
    /**
     * \brief Just a default constructor. My compiler made me do it!
     */
    Field() {
        name = "";
    };

    /**
     * \brief Creates a Field object with a given name
     * @param name The intended name for the Field
     */
    Field(std::string name) {
        setName(name);
    };

    /**
     * \brief Changes name of the Field
     * @param str New name for the Field
     */
    void setName(std::string str) {
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
    };

    /**
     * \brief Gives the data the Field holds
     * @return The data of the Field
     */
    T getData() {
        return data;
    };

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
        setData( std::stod(str) );
    };
};

/**
 * \brief Used to hold multiple real number values
 */
class DoubleVecField : public Field< std::vector<   double    > > {
public:
    using Field::Field;
    void parse(std::string str) {
        std::vector< double > vec;
        std::istringstream ss(str);
        std::string token;
        while ( getline(ss, token, ',') ) {
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
        str = (str.find_first_not_of(" \t\n"),str.find_last_not_of(" \t\n"));
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
    void registerSizeField(std::string name) {
        size_fields.push_back(SizeField(name));
        stringToFieldTypeMap.insert(std::pair<std::string, FieldType>(name, SIZE_FIELD));
    };

    /**
     * \brief Registers field into vector holding fields for multiple positive integer values
     * @param name
     */
    void registerSizeVecField(std::string name) {
        size_vec_fields.push_back(SizeVecField(name));
        stringToFieldTypeMap.insert(std::pair<std::string, FieldType>(name, SIZE_VEC_FIELD));
    };

    /**
     * \brief Registers field into vector holding fields for string values
     * @param name
     */
    void registerStringField(std::string name) {
        string_fields.push_back(StringField(name));
        stringToFieldTypeMap.insert(std::pair<std::string, FieldType>(name, STRING_FIELD));
    };

    /**
     * \brief Registers field into vector holding fields for multiple string values
     * @param name
     */
    void registerStringVecField(std::string name) {
        string_vec_fields.push_back(StringVecField(name));
        stringToFieldTypeMap.insert(std::pair<std::string, FieldType>(name, STRING_VEC_FIELD));
    };

    /**
     * \brief Registers field into vector holding fields for real number values
     * @param name
     */
    void registerDoubleField(std::string name) {
        double_fields.push_back(DoubleField(name));
        stringToFieldTypeMap.insert(std::pair<std::string, FieldType>(name, DOUBLE_FIELD));
    };

    /**
     * \brief Registers field into vector holding fields for multiple real number values
     * @param name
     */
    void registerDoubleVecField(std::string name) {
        double_vec_fields.push_back(DoubleVecField(name));
        stringToFieldTypeMap.insert(std::pair<std::string, FieldType>(name, DOUBLE_VEC_FIELD));
    };

    /**
     * Processes a full line of an input file after it has been stripped of all comments
     * @param line The line of the input file to process
     */
    void parseLine(std::string line);

    /**
     * \brief Prints all of the information pertaining to this Category
     */
    void printCategory();

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
    template <class T>
    typename std::vector<T>::iterator findNameInFields(std::string name, std::vector<T> &fields) {
        return std::find_if(fields.begin(), fields.end(),
                            [name](T & t) -> bool {
                                return t.getName().compare(name) == 0;
                            }
        );
    };

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

#endif //PNAB_FILEPARSER_H
