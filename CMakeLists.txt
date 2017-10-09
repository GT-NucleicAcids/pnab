cmake_minimum_required(VERSION 3.5)
project(pNAB)

set(CMAKE_CXX_STANDARD 11)

# set the module path
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)

# find and setup openbabel
find_package(OpenBabel2 REQUIRED)
include_directories(${OPENBABEL2_INCLUDE_DIR})

find_package(Eigen3 REQUIRED) # find and setup Eigen3 if available
set(EIGEN2_SUPPORT_STAGE10_FULL_EIGEN2_API TRUE)
include_directories(${EIGEN3_INCLUDE_DIR})

#find_package(OpenMP)
#if (OPENMP_FOUND)
#    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#endif()

set(SOURCE_FILES src/main.cpp src/Containers.h
        src/SearchTypes/SearchMethod.cpp src/SearchTypes/SearchMethod.h
                 src/UnitChain.cpp src/UnitChain.h src/ChainForceField.cpp
                 src/ChainForceField.h src/SearchTypes/RandomRotorSearch.cpp src/SearchTypes/RandomRotorSearch.h src/Containers.cpp src/FileParser.cpp src/FileParser.h)

add_executable(pNAB ${SOURCE_FILES})
target_link_libraries(pNAB ${OPENBABEL2_LIBRARIES})

# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxygen/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(doc
            ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating API documentation with Doxygen" VERBATIM
            )
endif(DOXYGEN_FOUND)