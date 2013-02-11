# NOTE: You will need to have the GNU environment loaded, e.g., via the command
#       
# module swap PrgEnv-intel PrgEnv-gnu
#

# The Cray wrappers
set(COMPILER_DIR /opt/cray/craype/1.01/bin)
set(CMAKE_C_COMPILER       ${COMPILER_DIR}/cc)
set(CMAKE_CXX_COMPILER     ${COMPILER_DIR}/CC)
set(CMAKE_Fortran_COMPILER ${COMPILER_DIR}/ftn)

# This is just a hack, as this machine always uses the above wrappers
set(MPI_C_COMPILER ${CMAKE_C_COMPILER})
set(MPI_CXX_COMPILER ${CMAKE_CXX_COMPILER})
set(MPI_Fortran_COMPILER ${CMAKE_Fortran_COMPILER})

set(C_FLAGS_PUREDEBUG "-g -static -Wl,-Bstatic")
set(C_FLAGS_PURERELEASE "-O3 -static -Wl,-Bstatic")
set(C_FLAGS_HYBRIDDEBUG "-g -static -Wl,-Bstatic")
set(C_FLAGS_HYBRIDRELEASE "-O3 -static -Wl,-Bstatic")
set(CXX_FLAGS_PUREDEBUG "-g -static -Wl,-Bstatic")
set(CXX_FLAGS_PURERELEASE "-O3 -static -Wl,-Bstatic")
set(CXX_FLAGS_HYBRIDDEBUG "-g -static -Wl,-Bstatic")
set(CXX_FLAGS_HYBRIDRELEASE "-O3 -static -Wl,-Bstatic")

set(OpenMP_CXX_FLAGS "-fopenmp")

#set(MATH_LIBS "/opt/xt-libsci/11.0.06/gnu/46/mc12/lib/libsci_gnu.a;/opt/gcc/4.7.1/snos/lib64/libgfortran.a;-lm")
#set(MATH_LIBS "/opt/cray/libsci/default/GNU/47/sandybridge/lib/libsci_gnu.a;/opt/gcc/4.7.2/snos/lib64/libgfortran.a;-lm")
#set(MATH_LIBS "/opt/cray/libsci/default/GNU/47/sandybridge/lib/libsci_gnu.a;/opt/gcc/4.7.2/snos/lib64/libgfortran.a;/opt/cray/xc-sysroot/5.0.15/usr/lib64/libm.a")
set(MATH_LIBS "/opt/cray/libsci/default/GNU/47/sandybridge/lib/libsci_gnu.a")
set(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
set(BUILD_SHARED_LIBS OFF)
set(CMAKE_EXE_LINKER_FLAGS "-static")
