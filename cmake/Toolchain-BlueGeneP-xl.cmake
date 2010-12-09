set(CMAKE_SYSTEM_NAME BlueGeneP-static)

# We need MPI C and CXX compilers, but only a serial Fortran compiler
set(CMAKE_C_COMPILER /bgsys/drivers/ppcfloor/comm/bin/mpixlc_r)
set(CMAKE_CXX_COMPILER /bgsys/drivers/ppcfloor/comm/bin/mpixlcxx_r)
set(CMAKE_Fortran_COMPILER /soft/apps/ibmcmp-aug2010/xlf/bg/11.1/bin/bgxlf_r) 

set(CXX_PURE_DEBUG_FLAGS "-g")
set(CXX_PURE_RELEASE_FLAGS "-g -O4")
set(CXX_HYBRID_DEBUG_FLAGS "-g")
set(CXX_HYBRID_RELEASE_FLAGS "-g -O4")

set(PTHREADS_C_FLAGS "-qthreaded")
set(OpenMP_CXX_FLAGS "-qsmp=omp:noauto -qthreaded")

# The remainder of the file is for linking BLAS/LAPACK functionality
set(ESSL_BASE "/soft/apps/ESSL-4.3.1-1")
set(IBMCMP_BASE "/soft/apps/ibmcmp-aug2010")
set(XLF_BASE "${IBMCMP_BASE}/xlf/bg/11.1/bglib")
set(XLSMP_BASE "${IBMCMP_BASE}/xlsmp/bg/1.7/bglib")

set(BGP_LAPACK "-L/soft/apps/LAPACK -llapack_bgp")
set(PURE_ESSL "-L${ESSL_BASE}/lib -lesslbg")
set(THREADED_ESSL "-L${ESSL_BASE}/lib -lesslsmpbg")
set(XLF_LIBS "-L${XLF_BASE} -lxlfmath -lxlf90_r")
set(XLOMP_SER "-L${XLSMP_BASE} -lxlomp_ser")
set(XLSMP "-L${XLSMP_BASE} -lxlsmp")

set(PURE_MATH_LIBS 
    "${BGP_LAPACK};${PURE_ESSL};${XLF_LIBS};${XLOMP_SER}")
set(HYBRID_MATH_LIBS 
    "${BGP_LAPACK};${THREADED_ESSL};${XLF_LIBS};${XLSMP}")

