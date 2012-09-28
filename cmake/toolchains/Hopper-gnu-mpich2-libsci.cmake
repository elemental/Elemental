# NOTE: You will need to have the GNU environment loaded, e.g., via the command
#       
#       module swap PrgEnv-pgi PrgEnv-gnu
#

# The Cray wrappers
set(COMPILER_DIR /opt/cray/xt-asyncpe/5.10/bin)
set(CMAKE_C_COMPILER       ${COMPILER_DIR}/cc)
set(CMAKE_CXX_COMPILER     ${COMPILER_DIR}/CC)
set(CMAKE_Fortran_COMPILER ${COMPILER_DIR}/ftn)

# This is just a hack, as this machine always uses the above wrappers
set(MPI_C_COMPILER ${CMAKE_C_COMPILER})
set(MPI_CXX_COMPILER ${CMAKE_CXX_COMPILER})
set(MPI_Fortran_COMPILER ${CMAKE_Fortran_COMPILER})

set(CXX_PURE_DEBUG_FLAGS "-g")
set(CXX_PURE_RELEASE_FLAGS "-O3 -ffast-math")
set(CXX_HYBRID_DEBUG_FLAGS "-g")
set(CXX_HYBRID_RELEASE_FLAGS "-O3 -ffast-math")

set(OpenMP_CXX_FLAGS "-fopenmp")

set(MATH_LIBS "/opt/xt-libsci/11.0.06/gnu/46/mc12/lib/libsci_gnu.a;/opt/gcc/4.7.1/snos/lib64/libgfortran.a;-lm")
