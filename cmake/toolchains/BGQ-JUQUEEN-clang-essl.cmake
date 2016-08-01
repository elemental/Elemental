# Set the system name so CMake uses the appropriate platform settings.
# NOTE: The platforms settings for BlueGeneP are the same for BlueGeneQ 
set(CMAKE_SYSTEM_NAME BlueGeneP-static)

set(GCC_ROOT  "/bgsys/local/gcc/4.8.1")
set(GCC_NAME  "powerpc64-bgq-linux")
set(CLANG_ROOT "/bgsys/local/clang/nightly")
set(CLANG_MPI_ROOT "${CLANG_ROOT}/mpi/bgclang")
set(IBMCMP_ROOT "/opt/ibmcmp")

set(BLAS_LIB "/bgsys/local/lib")
set(LAPACK_LIB "/bgsys/local/lapack/3.4.2/lib")
set(SCALAPACK_LIB "/bgsys/local/scalapack/2.0.2/lib")
set(ESSL_LIB "/bgsys/local/lib")

set(MPI_ROOT   "/bgsys/drivers/ppcfloor/comm/gcc")
set(PAMI_ROOT  "/bgsys/drivers/ppcfloor/comm/sys")
set(SPI_ROOT   "/bgsys/drivers/ppcfloor/spi")

# The serial compilers
set(CMAKE_C_COMPILER       "${CLANG_ROOT}/bin/bgclang")
set(CMAKE_CXX_COMPILER     "${CLANG_ROOT}/bin/bgclang++11")
set(CMAKE_Fortran_COMPILER "${GCC_ROOT}/bin/${GCC_NAME}-gfortran")

# The MPI wrappers for the C and C++ compilers
set(MPI_C_COMPILER   "${CLANG_MPI_ROOT}/bin/mpicc")
set(MPI_CXX_COMPILER "${CLANG_MPI_ROOT}/bin/mpic++11")

if(CMAKE_BUILD_TYPE MATCHES Debug)
  set(CXX_FLAGS "-g")
else()
  set(CXX_FLAGS "-O3")
endif()

set(CMAKE_THREAD_LIBS_INIT "-fopenmp")
set(OpenMP_CXX_FLAGS "-fopenmp")

set(XLF_LIB "${IBMCMP_ROOT}/xlf/bg/14.1/bglib64")
set(XLSMP_LIB "${IBMCMP_ROOT}/xlsmp/bg/3.1/bglib64")

if(EL_HYBRID)
  set(MATH_LIBS "-L${SCALAPACK_LIB} -lscalapack -L${ESSL_LIB} -lesslsmpbg -L${LAPACK_LIB} -llapack -L${ESSL_LIB} -lesslsmpbg -L${XLF_LIB} -lxlf90_r -L${XLSMP_LIB} -lxlsmp -lxlopt -lxlfmath -lxl -lpthread -ldl -Wl,--allow-multiple-definition")
else()
  set(MATH_LIBS "-L${SCALAPACK_LIB} -lscalapack -L${ESSL_LIB} -lesslbg -L${LAPACK_LIB} -llapack -L${ESSL_LIB} -lesslbg -L${XLF_LIB} -lxlf90_r -L${XLSMP_LIB} -lxlomp_ser -lxlopt -lxlfmath -lxl -lpthread -ldl -Wl,--allow-multiple-definition")
endif()
