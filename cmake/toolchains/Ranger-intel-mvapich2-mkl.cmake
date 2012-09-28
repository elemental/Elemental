# Set the serial Intel compilers
set(COMPILER_DIR /opt/apps/intel/10.1)
set(CMAKE_C_COMPILER       ${COMPILER_DIR}/cc/bin/icc)
set(CMAKE_CXX_COMPILER     ${COMPILER_DIR}/cc/bin/icpc)
set(CMAKE_Fortran_COMPILER ${COMPILER_DIR}/fc/bin/ifort)

# Set the MPI wrappers for the C and C++ compilers
set(MPI_COMPILER_DIR /opt/apps/intel10_1/mvapich2/1.2/bin)
set(MPI_C_COMPILER       ${MPI_COMPILER_DIR}/mpicc)
set(MPI_CXX_COMPILER     ${MPI_COMPILER_DIR}/mpicxx)
set(MPI_Fortran_COMPILER ${MPI_COMPILER_DIR}/mpif90)

set(CXX_PURE_DEBUG_FLAGS "-g")
set(CXX_PURE_RELEASE_FLAGS "-O3")
set(CXX_HYBRID_DEBUG_FLAGS "-g")
set(CXX_HYBRID_RELEASE_FLAGS "-O3")

set(OpenMP_CXX_FLAGS "-openmp")

set(MATH_LIBS 
    "-L/opt/apps/intel/mkl/10.0.1.014/lib/em64t -lmkl_em64t -lmkl -lguide -lpthread /opt/apps/intel/10.1/fc/lib/libifcore.a /opt/apps/intel/10.1/fc/lib/libsvml.a -lm")

