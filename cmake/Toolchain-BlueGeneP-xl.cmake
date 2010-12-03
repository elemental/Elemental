set(CMAKE_SYSTEM_NAME BlueGeneP-static)

# Serial compilers for Elemental+PMRRR
set(CMAKE_C_COMPILER /soft/apps/ibmcmp-aug2010/vac/bg/9.0/bin/bgxlc_r)
set(CMAKE_CXX_COMPILER /soft/apps/ibmcmp-aug2010/vacpp/bg/9.0/bin/bgxlC_r)
set(CMAKE_Fortran_COMPILER /soft/apps/ibmcmp-aug2010/xlf/bg/11.1/bin/bgxlf_r) 

# MPI compilers
set(MPI_CXX_COMPILER /bgsys/drivers/ppcfloor/comm/bin/mpixlcxx_r)
set(MPI_C_COMPILER /bgsys/drivers/ppcfloor/comm/bin/mpixlc_r)

set(CXX_PURE_DEBUG_FLAGS "-g")
set(CXX_PURE_RELEASE_FLAGS "-g -O4")
set(CXX_HYBRID_DEBUG_FLAGS "-g")
set(CXX_HYBRID_RELEASE_FLAGS "-g -O4")

set(MANUAL_OPENMP_CXX_FLAGS "-qsmp=omp:noauto -qthreaded")

# The remainder of the file is for linking BLAS/LAPACK functionality
set(ESSL_BASE "/soft/apps/ESSL-4.3.1-1")
set(IBMCMP_BASE "/soft/apps/ibmcmp-aug2010")
set(XLF_BASE "${IBMCMP_BASE}/xlf/bg/11.1/bglib")
set(XLSMP_BASE "${IBMCMP_BASE}/xlsmp/bg/1.7/bglib")

set(BGP_LAPACK "/soft/apps/LAPACK/liblapack_bgp.a")
set(PURE_ESSL "${ESSL_BASE}/lib/libesslbg.a")
set(THREADED_ESSL "${ESSL_BASE}/lib/libesslsmpbg.a")
set(XLF_LIBS "${XLF_BASE}/libxlfmath.a;${XLF_BASE}/libxlf90_r.a")
set(XLOMP_SER "${XLSMP_BASE}/libxlomp_ser.a")
set(XLSMP "${XLSMP_BASE}/libxlsmp.a")

set(PURE_BLAS_LAPACK_LIBS 
    "${BGP_LAPACK};${PURE_ESSL};${XLF_LIBS};${XLOMP_SER}")
set(HYBRID_BLAS_LAPACK_LIBS 
    "${BGP_LAPACK};${THREADED_ESSL};${XLF_LIBS};${XLSMP}")

