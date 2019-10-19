cmake_minimum_required(VERSION 3.5)
project(pnab)

set(CMAKE_CXX_STANDARD 11)
add_definitions(-DUSING_DYNAMIC_LIBS) # Necessary for openbabel in Windows

# set the module path
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)

find_package(OpenBabel3 REQUIRED)

find_package(Eigen3 REQUIRED) # find and setup Eigen3
set(EIGEN2_SUPPORT_STAGE10_FULL_EIGEN2_API TRUE)
include_directories(${EIGEN3_INCLUDE_DIR})

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

set(SOURCE_DIR pnab/src)
include_directories(${SOURCE_DIR})
set(SOURCE_FILES 
        ${SOURCE_DIR}/Containers.cpp
        ${SOURCE_DIR}/Chain.cpp
        ${SOURCE_DIR}/ConformationSearch.cpp
        )

# pybind11
find_package(pybind11 REQUIRED)
# Must add NO_EXTRAS; otherwise code does not work in Release mode
pybind11_add_module(bind NO_EXTRAS ${SOURCE_DIR}/binder.cpp ${SOURCE_FILES})
target_link_libraries(bind PRIVATE ${OPENBABEL3_LIBRARY})
target_include_directories(bind PRIVATE ${OPENBABEL3_INCLUDE_DIR})

add_library(restore_default_visibility INTERFACE)
target_compile_options(restore_default_visibility INTERFACE -fvisibility=default)
target_link_libraries(bind PRIVATE restore_default_visibility)
