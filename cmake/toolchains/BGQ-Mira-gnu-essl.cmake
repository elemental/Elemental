# Set the system name so CMake uses the appropriate platform settings.
# NOTE: The platforms settings for BlueGeneP are the same for BlueGeneQ 
set(CMAKE_SYSTEM_NAME BlueGeneP-static)

# The serial compilers are typically set here, but, at least at some point,
# it was easier to set them to the MPI compilers.
set(GCC_DIR    "/bgsys/drivers/ppcfloor/gnu-linux-4.7.2")
set(MPI_DIR    "/bgsys/drivers/ppcfloor/comm/gcc")

set(CMAKE_C_COMPILER       ${GCC_DIR}/bin/powerpc64-bgq-linux-gcc)
set(CMAKE_CXX_COMPILER     ${GCC_DIR}/bin/powerpc64-bgq-linux-g++)
set(CMAKE_Fortran_COMPILER ${GCC_DIR}/bin/powerpc64-bgq-linux-gfortran)

# The MPI wrappers for the C and C++ compilers
set(MPI_C_COMPILER   ${MPI_DIR}/bin/mpicc)
set(MPI_CXX_COMPILER ${MPI_DIR}/bin/mpicxx)

set(PAMI_ROOT "/bgsys/drivers/ppcfloor/comm/sys")
set(SPI_ROOT  "/bgsys/drivers/ppcfloor/spi")

if(CMAKE_BUILD_TYPE MATCHES Debug)
  set(CXX_FLAGS "-g")
else()
  set(CXX_FLAGS "-g -O3")
endif()

set(CMAKE_THREAD_LIBS_INIT "-fopenmp")
set(OpenMP_CXX_FLAGS "-fopenmp")

set(LAPACK_LIB "/soft/libraries/alcf/current/gcc/LAPACK/lib")
set(ESSL_LIB "/soft/libraries/essl/current/essl/5.1/lib64")
set(IBMCMP_ROOT "/soft/compilers/ibmcmp-nov2012")
set(XLF_LIB "${IBMCMP_ROOT}/xlf/bg/14.1/bglib64")
set(XLSMP_LIB "${IBMCMP_ROOT}/xlsmp/bg/3.1/bglib64")

set(LAPACK_FLAGS "-L${LAPACK_LIB} -llapack")
set(XLF_FLAGS "-L${XLF_LIB} -lxlf90_r")

if(EL_HYBRID)
  set(ESSL_FLAGS "-L${ESSL_LIB} -lesslsmpbg") 
  set(XL_FLAGS "-L${XLSMP_LIB} -lxlsmp")
else()
  set(ESSL_FLAGS "-L${ESSL_LIB} -lesslbg")
  set(XL_FLAGS "-L${XLSMP_LIB} -lxlomp_ser")
endif()

#set(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
#set(CMAKE_EXE_LINKER_FLAGS "-static")
set(MATH_LIBS "-L/soft/libraries/alcf/current/gcc/SCALAPACK/lib -lscalapack ${LAPACK_FLAGS} ${ESSL_FLAGS} ${XLF_FLAGS} ${XL_FLAGS} -lxlopt -lxlfmath -lxl -lgfortran -lm -lpthread -ldl -Wl,--allow-multiple-definition")
