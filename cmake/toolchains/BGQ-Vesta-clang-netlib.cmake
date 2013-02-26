#set(CMAKE_SYSTEM_NAME BlueGeneQ-static)

# The serial compilers
set(GCC_ROOT "/bgsys/drivers/ppcfloor/gnu-linux")
set(GCC_NAME "powerpc64-bgq-linux")
set(CLANG_ROOT "/home/projects/llvm")
set(CMAKE_C_COMPILER       "${CLANG_ROOT}/bin/bgclang")
set(CMAKE_CXX_COMPILER     "${CLANG_ROOT}/bin/bgclang++")
set(CMAKE_Fortran_COMPILER "${GCC_ROOT}/bin/${GCC_NAME}-gfortran")

# The MPI wrappers for the C and C++ compilers
set(MPI_ROOT  "/bgsys/drivers/ppcfloor/comm/gcc")
set(PAMI_ROOT "/bgsys/drivers/ppcfloor/comm/sys")
set(SPI_ROOT  "/bgsys/drivers/ppcfloor/spi")
set(MPI_C_COMPILE_FLAGS   "")
set(MPI_CXX_COMPILE_FLAGS "")
set(MPI_C_INCLUDE_PATH   "${MPI_ROOT}/include")
set(MPI_CXX_INCLUDE_PATH "${MPI_ROOT}/include")
set(MPI_C_LINK_FLAGS   "-L${MPI_ROOT}/lib -L${PAMI_ROOT}/lib -L${SPI_ROOT}/lib")
set(MPI_CXX_LINK_FLAGS ${MPI_C_LINK_FLAGS})
set(MPI_C_LIBRARIES "-lmpich -lopa -lmpl -ldl -lpami -lSPI -lSPI_cnk -lpthread -lrt -lstdc++")
set(MPI_CXX_LIBRARIES "-lcxxmpich ${MPI_C_LIBRARIES}")

#set(CXX_FLAGS_PUREDEBUG "-g -static -Bstatic")
#set(CXX_FLAGS_PURERELEASE "-g -O2 -static -Bstatic")
#set(CXX_FLAGS_HYBRIDDEBUG "-g -static -Bstatic")
#set(CXX_FLAGS_HYBRIDRELEASE "-g -O2 -static -Bstatic")
set(CXX_FLAGS_PUREDEBUG "-g")
set(CXX_FLAGS_PURERELEASE "-g -O2")
set(CXX_FLAGS_HYBRIDDEBUG "-g")
set(CXX_FLAGS_HYBRIDRELEASE "-g -O2")

#set(CMAKE_THREAD_LIBS_INIT "-fopenmp")
#set(OpenMP_CXX_FLAGS "-fopenmp")

##############################################################

# set the search path for the environment coming with the compiler
# and a directory where you can install your own compiled software
set(CMAKE_FIND_ROOT_PATH
    /bgsys/drivers/ppcfloor/
    /bgsys/drivers/ppcfloor/gnu-linux/powerpc64-bgq-linux
    /bgsys/drivers/ppcfloor/comm/gcc
    /bgsys/drivers/ppcfloor/comm/sys/
    /bgsys/drivers/ppcfloor/spi/)

# adjust the default behaviour of the FIND_XXX() commands:
# search headers and libraries in the target environment, search
# programs in the host environment
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)

##############################################################

#set(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
#set(CMAKE_EXE_LINKER_FLAGS "-static")
set(MATH_LIBS "-L/soft/libraries/alcf/current/gcc/LAPACK/lib -llapack -L/soft/libraries/alcf/current/gcc/BLAS/lib -lblas -lgfortran -lm")
