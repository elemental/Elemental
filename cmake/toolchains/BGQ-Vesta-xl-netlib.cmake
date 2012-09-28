set(CMAKE_SYSTEM_NAME BlueGeneQ-static)

set(IBMCMP_ROOT "/soft/compilers/ibmcmp-may2012")
set(MPI_ROOT    "/bgsys/drivers/ppcfloor/comm/xl")
set(PAMI_ROOT   "/bgsys/drivers/ppcfloor/comm/sys")
set(SPI_ROOT    "/bgsys/drivers/ppcfloor/spi")

# The serial XL compilers
#set(CMAKE_C_COMPILER       ${IBMCMP_ROOT}/vac/bg/12.1/bin/bgxlc_r)
#set(CMAKE_CXX_COMPILER     ${IBMCMP_ROOT}/vacpp/bg/12.1/bin/bgxlC_r)
#set(CMAKE_Fortran_COMPILER ${IBMCMP_ROOT}/xlf/bg/14.1/bin/bgxlf_r)
set(CMAKE_C_COMPILER       ${MPI_ROOT}/bin/mpixlc_r)
set(CMAKE_CXX_COMPILER     ${MPI_ROOT}/bin/mpixlcxx_r)
set(CMAKE_Fortran_COMPILER ${MPI_ROOT}/bin/mpixlf77_r)

# The MPI wrappers for the XL C and C++ compilers
set(MPI_C_COMPILER   ${MPI_ROOT}/bin/mpixlc_r)
set(MPI_CXX_COMPILER ${MPI_ROOT}/bin/mpixlcxx_r)

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

set(CXX_PURE_DEBUG_FLAGS "-g -O0 -qstrict -qnoipa")
set(CXX_PURE_RELEASE_FLAGS "-g -O3 -qarch=qp -qtune=qp -qsimd=auto
-qhot=level=1 -qprefetch -qunroll=yes -qreport -qnoipa")
set(CXX_HYBRID_DEBUG_FLAGS "-g -O0 -qstrict -qnoipa")
set(CXX_HYBRID_RELEASE_FLAGS "-g -O3 -qarch=qp -qtune=qp -qsimd=auto
-qhot=level=2 -qprefetch -qunroll=yes -qreport -qnoipa -qsmp=omp")

set(CMAKE_THREAD_LIBS_INIT "-qsmp=omp -qnoipa")
set(OpenMP_CXX_FLAGS "-qsmp=omp -qnoipa")

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

set(BLAS_LIB "/soft/libraries/alcf/current/xl/BLAS/lib")
set(LAPACK_LIB "/soft/libraries/alcf/current/xl/LAPACK/lib")
set(XLF_LIB "${IBMCMP_ROOT}/xlf/bg/14.1/bglib64")
set(XLSMP_LIB "${IBMCMP_ROOT}/xlsmp/bg/3.1/bglib64")
set(MATH_LIBS "-L${LAPACK_LIB} -llapack -L${BLAS_LIB} -lblas
-L${XLF_LIB} -lxlf90_r -L${XLSMP_LIB} -lxlsmp -lxlopt -lxlfmath -lxl
-lgfortran -lpthread -ldl")
