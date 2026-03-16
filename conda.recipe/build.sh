#!/bin/bash

set -e
set -x

BLD_DIR=`pwd`

SRC_DIR=$RECIPE_DIR/..

pushd $SRC_DIR
version=`git rev-parse --short HEAD`
popd

echo $version > __conda_version__.txt

export CMAKE_OSX_DEPLOYMENT_TARGET=""

cmake -D CMAKE_OSX_DEPLOYMENT_TARGET="" \
            -D CMAKE_INSTALL_PREFIX="${PREFIX}" \
            -D MPI_C_INCLUDE_PATH:STRING="${PREFIX}/include" \
            -D MPI_C_LINK_FLAGS:STRING="-L${PREFIX}/lib/" \
            -D MPI_C_LIBRARIES:STRING="-lmpi -lpmpi -L${PREFIX}/lib/" $SRC_DIR

patch -p0 < ${RECIPE_DIR}/cmake_parmetis.patch

#catch the failure here
make -j4 || :
cp download/parmetis/build/libparmetis/libparmetis.dylib ${PREFIX}/lib/libparmetis.dylib
cp download/parmetis/build/metis/libmetis/libmetis.dylib ${PREFIX}/lib/libmetis.dylib
make -j4

make install

cp -r ${PREFIX}/python/El ${PREFIX}/lib/python2.7/site-packages/
