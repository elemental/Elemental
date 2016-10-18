# The serial Intel compilers
set(COMPILER_DIR /opt/apps/intel/15/composer_xe_2015.2.164/bin/intel64)
set(CMAKE_C_COMPILER       ${COMPILER_DIR}/icc)
set(CMAKE_CXX_COMPILER     ${COMPILER_DIR}/icpc)
set(CMAKE_Fortran_COMPILER ${COMPILER_DIR}/ifort)

# The MPI wrappers for the C and C++ compilers
set(MPI_COMPILER_DIR /opt/apps/intel15/mvapich2/2.1/bin)
set(MPI_C_COMPILER       ${MPI_COMPILER_DIR}/mpicc)
set(MPI_CXX_COMPILER     ${MPI_COMPILER_DIR}/mpicxx)
set(MPI_Fortran_COMPILER ${MPI_COMPILER_DIR}/mpif90)

if(CMAKE_BUILD_TYPE MATCHES Debug)
  set(CXX_FLAGS "-g")
else()
  set(CXX_FLAGS "-O3")
endif()

set(OpenMP_CXX_FLAGS "-openmp")

set(MATH_LIBS "-L/opt/apps/intel/15/composer_xe_2015.2.164/mkl/lib/intel64 -L/opt/apps/intel/15/composer_xe_2015.2.164/compiler/lib/intel64 -lmkl_rt")
