#!/bin/bash

conda create -n pnab -c conda-forge python=3.7 numpy cmake openbabel eigen \
    pybind11 pyyaml py3dmol gcc_linux-64 gxx_linux-64 pytest pytest-cov codecov

source activate pnab

PREFIX=$CONDA_PREFIX
SRC_DIR=./
PY_ABBR=python3.7m
SP_DIR=$PREFIX/lib/python3.7/site-packages

${PREFIX}/bin/cmake \
    -H${SRC_DIR} \
    -Bbuild \
    -DCMAKE_INSTALL_PREFIX=${PREFIX} \
    -DCMAKE_BUILD_TYPE=Debug \
    -DCMAKE_C_COMPILER=${CC} \
    -DCMAKE_CXX_COMPILER=${CXX} \
    -DCMAKE_C_FLAGS="${CFLAGS}" \
    -DCMAKE_CXX_FLAGS="${CXXFLAGS}" \
    -DPYTHON_EXECUTABLE="${PREFIX}/bin/python" \
    -DPYTHON_LIBRARY="${PREFIX}/lib/lib${PY_ABBR}${SHLIB_EXT}" \
    -DPYTHON_INCLUDE_DIR="${PREFIX}/include/${PY_ABBR}" \
    -DBUILD_SHARED_LIBS=ON \
    -DBUILD_DOCS=OFF \
    -DENABLE_OPENMP=OFF \
    -DENABLE_XHOST=OFF \
    -DENABLE_GENERIC=ON \
    -Dopenbabel2_DIR="${PREFIX}/lib/cmake/openbabel2" \
    -Dpybind11_DIR="${PREFIX}/share/cmake/pybind11" \
    -DCMAKE_PREFIX_PATH="${PREFIX}"

# build
cd build
make -j${CPU_COUNT}
cp bind.*so ../pnab
cd ..

# install
cp -R pnab ${SP_DIR}
cp build/bind.*.so ${SP_DIR}/pnab
ls -l ${SP_DIR}/pnab

# test
pytest --cov
