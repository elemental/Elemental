# We need MPI C and CXX compilers, but only a serial Fortran compiler
set(CMAKE_C_COMPILER /opt/apps/gcc4_4/mvapich2/1.2/bin/mpicc)
set(CMAKE_CXX_COMPILER /opt/apps/gcc4_4/mvapich2/1.2/bin/mpicxx)
set(CMAKE_Fortran_COMPILER /opt/apps/gcc_amd/4.4.5/bin/gfortran)

set(CXX_PURE_DEBUG_FLAGS "-g")
set(CXX_PURE_RELEASE_FLAGS "-O3")
set(CXX_HYBRID_DEBUG_FLAGS "-g")
set(CXX_HYBRID_RELEASE_FLAGS "-O3")

set(OpenMP_CXX_FLAGS "-fopenmp")

set(MATH_LIBS 
    "-L/opt/apps/intel/mkl/10.0.1.014/lib/em64t -lmkl_em64t -lmkl -lguide -lpthread /opt/apps/gcc_amd/4.4.5/lib64/libgfortran.a -lm")

