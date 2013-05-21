# The serial Intel compilers
set(COMPILER_DIR /opt/apps/intel/13/composer_xe_2013.2.146/bin/intel64)
set(CMAKE_C_COMPILER       ${COMPILER_DIR}/icc)
set(CMAKE_CXX_COMPILER     ${COMPILER_DIR}/icpc)
set(CMAKE_Fortran_COMPILER ${COMPILER_DIR}/ifort)

# The MPI wrappers for the C and C++ compilers
set(MPI_COMPILER_DIR /opt/apps/intel13/mvapich2/1.9/bin)
set(MPI_C_COMPILER       ${MPI_COMPILER_DIR}/mpicc)
set(MPI_CXX_COMPILER     ${MPI_COMPILER_DIR}/mpicxx)
set(MPI_Fortran_COMPILER ${MPI_COMPILER_DIR}/mpif90)

set(MPI_C_COMPILE_FLAGS "")
set(MPI_CXX_COMPILE_FLAGS "")
set(MPI_Fortran_COMPILE_FLAGS "")
set(MPI_C_INCLUDE_PATH       /opt/apps/intel13/mvapich2/1.9/include)
set(MPI_CXX_INCLUDE_PATH     ${MPI_C_INCLUDE_PATH})
set(MPI_Fortran_INCLUDE_PATH ${MPI_C_INCLUDE_PATH})
set(MPI_C_LINK_FLAGS "-Wl,-rpath,/opt/apps/limic2/0.5.5/lib -L/opt/apps/limic2/0.5.5/lib -L/opt/apps/intel13/mvapich2/1.9/lib -L/opt/ofed/lib64/")
set(MPI_CXX_LINK_FLAGS ${MPI_C_LINK_FLAGS})
set(MPI_Fortran_LINK_FLAGS ${MPI_C_LINK_FLAGS})
set(MPI_BASE_LIBS 
    "-lmpich -lopa -llimic2 -lpthread -lrdmacm -libverbs -libumad -ldl -lrt")
set(MPI_C_LIBRARIES "-limf ${MPI_BASE_LIBS}")
set(MPI_CXX_LIBRARIES "-limf -lmpichcxx ${MPI_BASE_LIBS}")
set(MPI_Fortran_LIBRARIES "-limf -lmpichf90 ${MPI_BASE_LIBS}")

set(CXX_FLAGS_PUREDEBUG "-g")
set(CXX_FLAGS_PURERELEASE "-O3")
set(CXX_FLAGS_HYBRIDDEBUG "-g")
set(CXX_FLAGS_HYBRIDRELEASE "-O3")

set(OpenMP_CXX_FLAGS "-openmp")

set(MATH_LIBS "-mkl")
