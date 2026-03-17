#!/bin/sh

NPROC=`nproc`

git checkout tags/v0.87.5

mkdir build

cd build

cmake \
-DCMAKE_INSTALL_PREFIX="${PREFIX}" \
-DCMAKE_CXX_COMPILER="${PREFIX}/bin/g++" \
-DCMAKE_C_COMPILER="${PREFIX}/bin/gcc" \
-DCMAKE_Fortran_COMPILER="${PREFIX}/bin/gfortran" \
-DEL_USE_64BIT_INTS=ON \
-DEL_HAVE_QUADMATH=ON \
-DCMAKE_BUILD_TYPE=Release \
-DEL_HYBRID=OFF \
-DBUILD_SHARED_LIBS=ON \
-DMATH_LIBS="-L${PREFIX}/lib -lopenblas -lm" \
-DINSTALL_PYTHON_PACKAGE=ON \
-DEL_DISABLE_VALGRIND=ON \
-DEL_DISABLE_PARMETIS=ON \
.. 


make -j $NPROC

make install
