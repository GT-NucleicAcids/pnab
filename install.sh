#!/bin/bash

conda create -n pnab -c conda-forge python numpy cmake openbabel eigen \
    pybind11 pyyaml py3dmol gcc_linux-64 gxx_linux-64 pytest graphviz

source activate pnab

cmake -Bbuild -DCMAKE_BUILD_TYPE=Release -DBUILD_DOCS=OFF

# build
cd build
make
cp bind.*so ../pnab
cd ..

# install
SP_DIR=$CONDA_PREFIX/lib/python3.7/site-packages
cp -R pnab ${SP_DIR}
cp -R tests ${SP_DIR}/pnab
cp build/bind.*.so ${SP_DIR}/pnab
ls -l ${SP_DIR}/pnab

# test
pytest
