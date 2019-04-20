
if [ "$(uname)" == "Linux" ]; then

    # configure
    ${BUILD_PREFIX}/bin/cmake \
        -H${SRC_DIR} \
        -Bbuild3 \
        -DCMAKE_INSTALL_PREFIX=${PREFIX} \
        -DCMAKE_BUILD_TYPE=Debug \
        -DCMAKE_C_COMPILER=${CC} \
        -DCMAKE_CXX_COMPILER=${CXX} \
        -DCMAKE_C_FLAGS="${CFLAGS}" \
        -DCMAKE_CXX_FLAGS="${CXXFLAGS}" \
        -DPYTHON_EXECUTABLE="${PYTHON}" \
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

fi

# build
cd build3
make -j${CPU_COUNT}
cd ..

# install
cp -R pNAB ${SP_DIR}
cp build3/bind.*.so ${SP_DIR}/pNAB
ls -l ${SP_DIR}/pNAB
