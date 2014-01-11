#set(CMAKE_SYSTEM_NAME BlueGeneQ-static)

set(GCC_ROOT  "/bgsys/drivers/ppcfloor/gnu-linux")
set(GCC_NAME  "powerpc64-bgq-linux")
set(CLANG_ROOT "/home/projects/llvm")
set(CLANG_MPI_ROOT "/home/projects/llvm/mpi/bgclang")
set(IBMCMP_ROOT "$ENV{IBM_MAIN_DIR}")

set(BLAS_LIB "/soft/libraries/alcf/current/xl/BLAS/lib")
set(LAPACK_LIB "/soft/libraries/alcf/current/xl/LAPACK/lib")
set(ESSL_LIB "/soft/libraries/essl/current/essl/5.1/lib64")

set(MPI_ROOT   "/bgsys/drivers/ppcfloor/comm/gcc")
set(PAMI_ROOT  "/bgsys/drivers/ppcfloor/comm/sys")
set(SPI_ROOT   "/bgsys/drivers/ppcfloor/spi")

# The serial compilers
set(CMAKE_C_COMPILER       "${CLANG_ROOT}/wbin/bgclang")
set(CMAKE_CXX_COMPILER     "${CLANG_ROOT}/wbin/bgclang++11")
set(CMAKE_Fortran_COMPILER "${GCC_ROOT}/bin/${GCC_NAME}-gfortran")

# The MPI wrappers for the C and C++ compilers
set(MPI_C_COMPILER   "${CLANG_MPI_ROOT}/bin/mpiclang")
set(MPI_CXX_COMPILER "${CLANG_MPI_ROOT}/bin/mpiclang++11")

set(MPI_C_COMPILE_FLAGS    "")
set(MPI_CXX_COMPILE_FLAGS  "")
set(MPI_C_INCLUDE_PATH     "${MPI_ROOT}/include")
set(MPI_CXX_INCLUDE_PATH   "${MPI_ROOT}/include")
set(MPI_C_LINK_FLAGS       "-L${MPI_ROOT}/lib -L${PAMI_ROOT}/lib -L${SPI_ROOT}/lib")
set(MPI_CXX_LINK_FLAGS     "${MPI_C_LINK_FLAGS}")
set(MPI_C_LIBRARIES       "${MPI_C_LINK_FLAGS}   -lmpich -lopa -lmpl -lpami -lSPI -lSPI_cnk -lrt -lpthread -lstdc++ -lpthread")
set(MPI_CXX_LIBRARIES     "${MPI_CXX_LINK_FLAGS} -lcxxmpich ${MPI_C_LIBRARIES}")

if(CMAKE_BUILD_TYPE MATCHES PureDebug OR
   CMAKE_BUILD_TYPE MATCHES HybridDebug)
  set(CXX_FLAGS "-g")
else()
  set(CXX_FLAGS "-O3")
endif()

set(CMAKE_THREAD_LIBS_INIT "-fopenmp")
set(OpenMP_CXX_FLAGS "-fopenmp")

##############################################################

# set the search path for the environment coming with the compiler
# and a directory where you can install your own compiled software
set(CMAKE_FIND_ROOT_PATH
    /bgsys/drivers/ppcfloor
    /bgsys/drivers/ppcfloor/gnu-linux/powerpc64-bgq-linux
    /bgsys/drivers/ppcfloor/comm/xl
    /bgsys/drivers/ppcfloor/comm/sys
    /bgsys/drivers/ppcfloor/spi
)

# adjust the default behaviour of the FIND_XXX() commands:
# search headers and libraries in the target environment, search
# programs in the host environment
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)

##############################################################

set(XLF_LIB "${IBMCMP_ROOT}/xlf/bg/14.1/bglib64")
set(XLSMP_LIB "${IBMCMP_ROOT}/xlsmp/bg/3.1/bglib64")
if(CMAKE_BUILD_TYPE MATCHES PureDebug OR
   CMAKE_BUILD_TYPE MATCHES PureRelease)
  set(MATH_LIBS "-L${ESSL_LIB} -lesslbg -L${LAPACK_LIB} -llapack -L${ESSL_LIB} -lesslbg -L${XLF_LIB} -lxlf90_r -L${XLSMP_LIB} -lxlomp_ser -lxlopt -lxlfmath -lxl -lpthread -ldl -Wl,--allow-multiple-definition")
else()
  set(MATH_LIBS "-L${ESSL_LIB} -lesslsmpbg -L${LAPACK_LIB} -llapack -L${ESSL_LIB} -lesslsmpbg -L${XLF_LIB} -lxlf90_r -L${XLSMP_LIB} -lxlsmp -lxlopt -lxlfmath -lxl -lpthread -ldl -Wl,--allow-multiple-definition")
endif()
