# The serial GNU compilers
set(COMPILER_DIR /share/sw/software/GCC/4.8.2/bin)
set(CMAKE_C_COMPILER       ${COMPILER_DIR}/gcc)
set(CMAKE_CXX_COMPILER     ${COMPILER_DIR}/g++)
set(CMAKE_Fortran_COMPILER ${COMPILER_DIR}/gfortran)

# The MPI wrappers for the C and C++ compilers
set(MPI_COMPILER_DIR /share/sw/software/OpenMPI/1.7.3-GCC-4.8.2/bin)
set(MPI_C_COMPILER       ${MPI_COMPILER_DIR}/mpicc)
set(MPI_CXX_COMPILER     ${MPI_COMPILER_DIR}/mpicxx)
set(MPI_Fortran_COMPILER ${MPI_COMPILER_DIR}/mpif90)

if(CMAKE_BUILD_TYPE MATCHES Debug)
  set(CXX_FLAGS "-g")
else()
  set(CXX_FLAGS "-O3")
endif()

set(OpenMP_CXX_FLAGS "-fopenmp")

set(MATH_LIBS "-L/share/sw/software/imkl/11.1.2.144/composer_xe_2013_sp1.2.144/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_lapack95_lp64 -lmkl_blas95_lp64 -lmkl_sequential -lpthread /share/sw/software/GCC/4.8.2/lib64/libgfortran.a /share/sw/software/imkl/11.1.2.144/composer_xe_2013_sp1.2.144/mkl/lib/intel64/libimf.a /share/sw/software/imkl/11.1.2.144/composer_xe_2013_sp1.2.144/mkl/lib/intel64/libirc.a -lm")
