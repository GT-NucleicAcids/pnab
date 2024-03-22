#!/bin/bash

conda create -n pnab -c conda-forge python=3.10 numpy cmake openbabel eigen \
    "pybind11<2.10" pyyaml "nglview>=2.7, <3" gcc_linux-64 gxx_linux-64 pytest graphviz "ipywidgets>=7.5, <8" jupyterlab=3

source $CONDA_PREFIX/bin/activate pnab

cmake -Bbuild -DCMAKE_BUILD_TYPE=Release -DBUILD_DOCS=ON -DPYTHON_EXECUTABLE="$CONDA_PREFIX/bin/python" -DOPENBABEL_DIR="$CONDA_PREFIX"

# build
cd build
make
cp bind.*so ../pnab
cd ..

# install
SP_DIR=$CONDA_PREFIX/lib/python3.10/site-packages
cp -R pnab ${SP_DIR}
cp -R tests ${SP_DIR}/pnab
cp build/bind.*.so ${SP_DIR}/pnab
ls -l ${SP_DIR}/pnab

# test
pytest -s
