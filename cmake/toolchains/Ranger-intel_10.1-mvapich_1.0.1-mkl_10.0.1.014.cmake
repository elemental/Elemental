# We need MPI C and CXX compilers, but only a serial Fortran compiler
set(CMAKE_C_COMPILER /opt/apps/intel10_1/mvapich/1.0.1/bin/mpicc)
set(CMAKE_CXX_COMPILER /opt/apps/intel10_1/mvapich/1.0.1/bin/mpicxx)
set(CMAKE_Fortran_COMPILER /opt/apps/intel/10.1/fc/bin/ifort)

set(CXX_PURE_DEBUG_FLAGS "-g")
set(CXX_PURE_RELEASE_FLAGS "-O3")
set(CXX_HYBRID_DEBUG_FLAGS "-g")
set(CXX_HYBRID_RELEASE_FLAGS "-O3")

set(OpenMP_CXX_FLAGS "-openmp")

set(MATH_LIBS 
    "-L/opt/apps/intel/mkl/10.0.1.014/lib/em64t -lmkl_em64t -lmkl -lguide -lpthread /opt/apps/intel/10.1/fc/lib/libifcore.a /opt/apps/intel/10.1/fc/lib/libsvml.a -lm")

