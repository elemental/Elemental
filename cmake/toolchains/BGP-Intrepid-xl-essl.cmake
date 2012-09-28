set(CMAKE_SYSTEM_NAME BlueGeneP-static)

# The serial XL compilers
set(BGXLC_BASE /soft/apps/ibmcmp-aug2011/vacpp/bg/9.0)
set(BGXLF_BASE /soft/apps/ibmcmp-aug2011/xlf/bg/11.1)
set(CMAKE_C_COMPILER       ${BGXLC_BASE}/bin/bgxlc_r)
set(CMAKE_CXX_COMPILER     ${BGXLC_BASE}/bin/bgxlC_r)
set(CMAKE_Fortran_COMPILER ${BGXLF_BASE}/bin/bgxlf_r)

# The MPI wrappers for the XL C and C++ compilers
set(BGMPI_BASE /bgsys/drivers/ppcfloor/comm)
set(MPI_C_COMPILER   ${BGMPI_BASE}/bin/mpixlc_r)
set(MPI_CXX_COMPILER ${BGMPI_BASE}/bin/mpixlcxx_r)
set(MPI_Fortran_COMPILER ${BGMPI_BASE}/bin/mpixlf90_r)

set(CXX_PURE_DEBUG_FLAGS "-g")
set(CXX_PURE_RELEASE_FLAGS "-g -O3")
set(CXX_HYBRID_DEBUG_FLAGS "-g")
set(CXX_HYBRID_RELEASE_FLAGS "-g -O3")

set(CMAKE_THREAD_LIBS_INIT "-qthreaded")
set(OpenMP_CXX_FLAGS "-qsmp=omp:noauto -qthreaded")

# The remainder of the file is for linking BLAS/LAPACK functionality
set(ESSL_BASE "/soft/apps/ESSL-4.3.1-1")
set(IBMCMP_BASE "/soft/apps/ibmcmp-aug2011")
set(XLF_BASE "${IBMCMP_BASE}/xlf/bg/11.1/bglib")
set(XLSMP_BASE "${IBMCMP_BASE}/xlsmp/bg/1.7/bglib")

set(BGP_LAPACK "-L/soft/apps/LAPACK -llapack_bgp")
set(PURE_ESSL "-L${ESSL_BASE}/lib -lesslbg")
set(THREADED_ESSL "-L${ESSL_BASE}/lib -lesslsmpbg")
set(XLF_LIBS "-L${XLF_BASE} -lxlfmath -lxlf90_r")
set(XLOMP_SER "-L${XLSMP_BASE} -lxlomp_ser")
set(XLSMP "-L${XLSMP_BASE} -lxlsmp")

if(CMAKE_BUILD_TYPE MATCHES PureDebug OR
   CMAKE_BUILD_TYPE MATCHES PureRelease)
  set(MATH_LIBS "${BGP_LAPACK};${PURE_ESSL};${XLF_LIBS};${XLOMP_SER}")
else()
  set(MATH_LIBS "${BGP_LAPACK};${THREADED_ESSL};${XLF_LIBS};${XLSMP}")
endif()
