# The serial XL compilers
set(CMAKE_C_COMPILER /usr/bin/bgxlc_r)
set(CMAKE_CXX_COMPILER /usr/bin/bgxlc++_r)
set(CMAKE_Fortran_COMPILER /usr/bin/bgxlf_r)

# The MPI wrappers for the XL compilers
set(MPI_ROOT /bgsys/drivers/ppcfloor/comm/xl)
set(MPI_C_COMPILER ${MPI_ROOT}/bin/mpixlc_r)
set(MPI_CXX_COMPILER ${MPI_ROOT}/bin/mpixlcxx_r)
set(MPI_Fortran_COMPILER ${MPI_ROOT}/bin/mpixlf90_r)

set(MPI_C_COMPILE_FLAGS       "")
set(MPI_CXX_COMPILE_FLAGS     "")
set(MPI_Fortran_COMPILE_FLAGS "")
set(MPI_C_INCLUDE_PATH       "/bgsys/drivers/V1R1M1/ppc64/comm/sys/include;/bgsys/drivers/V1R1M1/ppc64;/bgsys/drivers/V1R1M1/ppc64/spi/include;/bgsys/drivers/V1R1M1/ppc64/spi/include/kernel/cnk;/bgsys/drivers/V1R1M1/ppc64/comm/xl/include")
set(MPI_CXX_INCLUDE_PATH     ${MPI_C_INCLUDE_PATH})
set(MPI_Fortran_INCLUDE_PATH ${MPI_C_INCLUDE_PATH})
set(MPI_C_LINK_FLAGS "-L/bgsys/drivers/V1R1M1/ppc64/comm/xl/lib -L/bgsys/drivers/V1R1M1/ppc64/comm/sys/lib -L/bgsys/drivers/V1R1M1/ppc64/spi/lib")
set(MPI_CXX_LINK_FLAGS ${MPI_C_LINK_FLAGS})
set(MPI_Fortran_LINK_FLAGS ${MPI_C_LINK_FLAGS})
set(MPI_BASE_LIBS "-lmpich -lopa -lmpl -lrt -ldl -lpami -lSPI -lSPI_cnk -lpthread -lrt -lstdc++")
set(MPI_C_LIBRARIES ${MPI_BASE_LIBS})
set(MPI_CXX_LIBRARIES "-lcxxmpich ${MPI_BASE_LIBS}")
set(MPI_Fortran_LIBRARIES ${MPI_BASE_LIBS})

set(CXX_FLAGS_PUREDEBUG "-g")
set(CXX_FLAGS_PURERELEASE "-g -O2")
set(CXX_FLAGS_HYBRIDDEBUG "-g")
set(CXX_FLAGS_HYBRIDRELEASE "-g -O2")

set(CMAKE_THREAD_LIBS_INIT "-qthreaded")
set(OpenMP_CXX_FLAGS "-qsmp=omp:noauto -qthreaded")

##############################################################

# set the search path for the environment coming with the compiler
# and a directory where you can install your own compiled software
set(CMAKE_FIND_ROOT_PATH
    /bgsys/drivers/ppcfloor/
    /bgsys/drivers/ppcfloor/comm/xl
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

if(CMAKE_BUILD_TYPE MATCHES PureDebug OR
   CMAKE_BUILD_TYPE MATCHES PureRelease)
  set(MATH_LIBS "-L/bgsys/local/lapack/3.3.0/lib -llapack -L/opt/ibmmath/essl/5.1/lib64 -lesslbg -L/opt/ibmcmp/xlf/bg/14.1/bglib64 -lxlfmath -lxlf90_r -lm")
else()
  set(MATH_LIBS "-L/bgsys/local/lapack/3.3.0/lib -llapack -L/opt/ibmmath/essl/5.1/lib64 -lesslsmpbg -L/opt/ibmcmp/xlf/bg/14.1/bglib64 -lxlfmath -lxlf90_r -lm")
endif()
