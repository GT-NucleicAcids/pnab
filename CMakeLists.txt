cmake_minimum_required(VERSION 3.5)
project(pNAB)

set(CMAKE_CXX_STANDARD 14)

# set the module path
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)

# find and setup openbabel
find_package(OpenBabel2 REQUIRED)
include_directories(${OPENBABEL2_INCLUDE_DIR})

find_package(Eigen3 REQUIRED) # find and setup Eigen3 if available
set(EIGEN2_SUPPORT_STAGE10_FULL_EIGEN2_API TRUE)
include_directories(${EIGEN3_INCLUDE_DIR})

set(SOURCE_FILES src/main.cpp src/Containers.cpp src/Containers.h)
add_executable(pNAB ${SOURCE_FILES})