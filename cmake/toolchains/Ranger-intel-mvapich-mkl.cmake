# Set the serial Intel compilers
set(COMPILER_DIR /opt/apps/intel/10.1)
set(CMAKE_C_COMPILER       ${COMPILER_DIR}/cc/bin/icc)
set(CMAKE_CXX_COMPILER     ${COMPILER_DIR}/cc/bin/icpc)
set(CMAKE_Fortran_COMPILER ${COMPILER_DIR}/fc/bin/ifort)

# Set the MPI wrappers for the C and C++ compilers
set(MPI_COMPILER_DIR /opt/apps/intel10_1/mvapich/1.0.1/bin)
set(MPI_C_COMPILER       ${MPI_COMPILER_DIR}/mpicc)
set(MPI_CXX_COMPILER     ${MPI_COMPILER_DIR}/mpicxx)
set(MPI_Fortran_COMPILER ${MPI_COMPILER_DIR}/mpif90)

set(MPI_C_COMPILE_FLAGS "")
set(MPI_CXX_COMPILE_FLAGS "")
set(MPI_Fortran_COMPILE_FLAGS "")
set(MPI_C_INCLUDE_PATH   /opt/apps/intel10_1/mvapich/1.0.1/include)
set(MPI_CXX_INCLUDE_PATH ${MPI_C_INCLUDE_PATH})
set(MPI_Fortran_INCLUDE_PATH ${MPI_C_INCLUDE_PATH})
set(MPI_C_LINK_FLAGS "-Wl,-rpath,/opt/apps/intel10_1/mvapich/1.0.1/lib/shared -Wl,-rpath,/opt/apps/intel10_1/mvapich/1.0.1/lib -i-dynamic -Wl,-rpath,/opt/apps/intel/10.1/fc/lib -Wl,-rpath,/opt/apps/intel/10.1/cc/lib -i-dynamic -L/opt/apps/intel10_1/mvapich/1.0.1/lib -L/opt/apps/intel10_1/mvapich/1.0.1/lib/shared  -L/usr/lib64 -L/opt/ofed/lib64")
set(MPI_CXX_LINK_FLAGS ${MPI_C_LINK_FLAGS})
set(MPI_Fortran_LINK_FLAGS ${MPI_C_LINK_FLAGS})
set(MPI_BASE_LIBS "-lmpich -libverbs -libumad -lpthread -lrt")
set(MPI_C_LIBRARIES "${MPI_BASE_LIBS}")
set(MPI_CXX_LIBRARIES "-lpmpich++ ${MPI_BASE_LIBS}")
set(MPI_Fortran_LIBRARIES "-lmpichf90nc -lmpichfarg ${MPI_BASE_LIBS}")

set(CXX_PUREDEBUG_FLAGS "-g")
set(CXX_PURERELEASE_FLAGS "-O3")
set(CXX_HYBRIDDEBUG_FLAGS "-g")
set(CXX_HYBRIDRELEASE_FLAGS "-O3")

set(OpenMP_CXX_FLAGS "-openmp")

set(MATH_LIBS 
    "-L/opt/apps/intel/mkl/10.0.1.014/lib/em64t -lmkl_em64t -lmkl -lguide -lpthread /opt/apps/intel/10.1/fc/lib/libifcore.a /opt/apps/intel/10.1/fc/lib/libsvml.a -lm")
