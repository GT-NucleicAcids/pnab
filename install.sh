#!/bin/bash

conda create -n pnab -c conda-forge python numpy cmake pkg-config openbabel eigen \
    pybind11 pyyaml nglview gcc_linux-64 gxx_linux-64 pytest graphviz

source activate pnab

cmake -Bbuild -DCMAKE_BUILD_TYPE=Release \
              -DBUILD_DOCS=ON \
              -DOPENBABEL3_INCLUDE_DIR=$CONDA_PREFIX/include/openbabel3 \
              -DOPENBABEL3_LIBRARY=$CONDA_PREFIX/lib/libopenbabel.so \
              -DOPENBABEL3_VERSION_MET=TRUE \
              -DPYTHON_EXECUTABLE="$CONDA_PREFIX/bin/python"

# build
cd build
make
cp bind.*so ../pnab
cd ..

# install
SP_DIR=$CONDA_PREFIX/lib/python3.8/site-packages
cp -R pnab ${SP_DIR}
cp -R tests ${SP_DIR}/pnab
cp build/bind.*.so ${SP_DIR}/pnab
ls -l ${SP_DIR}/pnab

# test
pytest
