set(CMAKE_SYSTEM_NAME BlueGeneQ-static)

set(XLF_ROOT  "/opt/ibmcmp/xlf/bg/14.1")
set(GCC_ROOT  "/bgsys/drivers/ppcfloor/gnu-linux")
set(GCC_NAME  "powerpc64-bgq-linux")
set(MPI_ROOT  "/bgsys/drivers/ppcfloor/comm/gcc")
set(PAMI_ROOT "/bgsys/drivers/ppcfloor/comm/sys")
set(SPI_ROOT  "/bgsys/drivers/ppcfloor/spi")

# The serial compilers
#set(CMAKE_C_COMPILER       ${GCC_ROOT}/bin/${GCC_NAME}-gcc)
#set(CMAKE_CXX_COMPILER     ${GCC_ROOT}/bin/${GCC_NAME}-g++)
#set(CMAKE_Fortran_COMPILER ${GCC_ROOT}/bin/${GCC_NAME}-gfortran)
set(CMAKE_C_COMPILER        ${MPI_ROOT}/bin/mpicc)
set(CMAKE_CXX_COMPILER      ${MPI_ROOT}/bin/mpicxx)
set(CMAKE_Fortran_COMPILER  ${MPI_ROOT}/bin/mpif77)

# The MPI wrappers for the C and C++ compilers
set(MPI_C_COMPILER   ${MPI_ROOT}/bin/mpicc)
set(MPI_CXX_COMPILER ${MPI_ROOT}/bin/mpicxx)

set(MPI_C_COMPILE_FLAGS   "")
set(MPI_CXX_COMPILE_FLAGS "")
set(MPI_C_INCLUDE_PATH   "${MPI_ROOT}/include")
set(MPI_CXX_INCLUDE_PATH "${MPI_ROOT}/include")
set(MPI_C_LINK_FLAGS   "-L${MPI_ROOT}/lib -L${PAMI_ROOT}/lib -L${SPI_ROOT}/lib")
set(MPI_CXX_LINK_FLAGS "-L${MPI_ROOT}/lib -L${PAMI_ROOT}/lib -L${SPI_ROOT}/lib")
set(MPI_C_LIBRARIES              "-lmpich -lopa -lmpl -lrt -ldl -lpami
-lSPI -lSPI_cnk -lpthread -lrt -lstdc++")
set(MPI_CXX_LIBRARIES "-lcxxmpich -lmpich -lopa -lmpl -lrt -ldl -lpami
-lSPI -lSPI_cnk -lpthread -lrt -lstdc++")

set(CXX_PURE_DEBUG_FLAGS "-g")
set(CXX_PURE_RELEASE_FLAGS "-g -O2")
set(CXX_HYBRID_DEBUG_FLAGS "-g")
set(CXX_HYBRID_RELEASE_FLAGS "-g -O2")

set(CMAKE_THREAD_LIBS_INIT "-fopenmp")
set(OpenMP_CXX_FLAGS "-fopenmp")

##############################################################

# set the search path for the environment coming with the compiler
# and a directory where you can install your own compiled software
set(CMAKE_FIND_ROOT_PATH
    /bgsys/drivers/ppcfloor/
    /bgsys/drivers/ppcfloor/gnu-linux/powerpc64-bgq-linux
    /bgsys/drivers/ppcfloor/comm/gcc
    /bgsys/drivers/ppcfloor/comm/sys/
    /bgsys/drivers/ppcfloor/spi/
)

# adjust the default behaviour of the FIND_XXX() commands:
# search headers and libraries in the target environment, search
# programs in the host environment
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)

##############################################################

set(MATH_LIBS "-L/soft/libraries/alcf/current/gcc/LAPACK/lib -llapack
-L/soft/libraries/alcf/current/gcc/BLAS/lib -lblas -lgfortran -lm")
