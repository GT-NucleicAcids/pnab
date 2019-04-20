cmake_minimum_required(VERSION 3.5)
project(pNAB)

set(CMAKE_CXX_STANDARD 11)

# set the module path
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)

# find and setup openbabel
find_package(OpenBabel2 CONFIG REQUIRED)

find_package(Eigen3 REQUIRED) # find and setup Eigen3 if available
set(EIGEN2_SUPPORT_STAGE10_FULL_EIGEN2_API TRUE)
include_directories(${EIGEN3_INCLUDE_DIR})

#find_package(OpenMP)
#if (OPENMP_FOUND)
#    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#endif()

# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(BUILD_DOCS AND DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxygen/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(doc ALL
            COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating API documentation with Doxygen" VERBATIM
            )
endif()

set(SOURCE_DIR pNAB/src)
include_directories(${SOURCE_DIR})
set(SOURCE_FILES 
        ${SOURCE_DIR}/Containers.cpp
        ${SOURCE_DIR}/Chain.cpp
        ${SOURCE_DIR}/SearchTypes/MonteCarloRotorSearch.cpp
        )

# pybind11
find_package(pybind11 REQUIRED)
pybind11_add_module(bind ${SOURCE_DIR}/binder.cpp ${SOURCE_FILES})
target_link_libraries(bind PRIVATE openbabel)
target_include_directories(bind PRIVATE ${OpenBabel2_INCLUDE_DIRS})

add_library(restore_default_visibility INTERFACE)
target_compile_options(restore_default_visibility INTERFACE -fvisibility=default)
target_link_libraries(bind PRIVATE restore_default_visibility)
