# The serial Intel compilers
set(COMPILER_DIR /opt/apps/intel/11.1/bin/intel64)
set(CMAKE_C_COMPILER       ${COMPILER_DIR}/icc)
set(CMAKE_CXX_COMPILER     ${COMPILER_DIR}/icpc)
set(CMAKE_Fortran_COMPILER ${COMPILER_DIR}/ifort)

# The MPI wrappers for the C and C++ compilers
set(MPI_COMPILER_DIR /opt/apps/intel11_1/mvapich2/1.6/bin)
set(MPI_C_COMPILER       ${MPI_COMPILER_DIR}/mpicc)
set(MPI_CXX_COMPILER     ${MPI_COMPILER_DIR}/mpicxx)
set(MPI_Fortran_COMPILER ${MPI_COMPILER_DIR}/mpif90)

set(MPI_C_COMPILE_FLAGS "")
set(MPI_CXX_COMPILE_FLAGS "")
set(MPI_Fortran_COMPILE_FLAGS "")
set(MPI_C_INCLUDE_PATH       /opt/apps/intel11_1/mvapich2/1.6/include)
set(MPI_CXX_INCLUDE_PATH     ${MPI_C_INCLUDE_PATH})
set(MPI_Fortran_INCLUDE_PATH ${MPI_C_INCLUDE_PATH})
set(MPI_C_LINK_FLAGS "-Wl,-rpath,/opt/apps/intel/11.1/lib/intel64 -L/opt/apps/intel/11.1/lib/intel64 -Wl,-rpath,/opt/apps/intel/11.1/lib/intel64 -Wl,-rpath,/opt/apps/limic2/0.5.4/lib -L/opt/apps/limic2/0.5.4/lib -L/opt/apps/intel11_1/mvapich2/1.6/lib -L/opt/apps/limic2/0.5.4/lib -L/opt/ofed/lib64/")
set(MPI_CXX_LINK_FLAGS ${MPI_C_LINK_FLAGS})
set(MPI_Fortran_LINK_FLAGS ${MPI_C_LINK_FLAGS})
set(MPI_BASE_LIBS 
    "-lmpich -lopa -llimic2 -lpthread -lrdmacm -libverbs -libumad -ldl -lrt")
set(MPI_C_LIBRARIES "-limf ${MPI_BASE_LIBS}")
set(MPI_CXX_LIBRARIES "-limf -lmpichcxx ${MPI_BASE_LIBS}")
set(MPI_Fortran_LIBRARIES "-limf -lmpichf90 ${MPI_BASE_LIBS}")

set(CXX_PURE_DEBUG_FLAGS "-g")
set(CXX_PURE_RELEASE_FLAGS "-O3")
set(CXX_HYBRID_DEBUG_FLAGS "-g")
set(CXX_HYBRID_RELEASE_FLAGS "-O3")

set(OpenMP_CXX_FLAGS "-openmp")

if(CMAKE_BUILD_TYPE MATCHES PureDebug OR
   CMAKE_BUILD_TYPE MATCHES PureRelease)
  set(MATH_LIBS "-L/opt/apps/intel/11.1/mkl/lib/em64t -lmkl_intel_lp64 -lmkl_sequential -lmkl_core /opt/apps/intel/11.1/lib/intel64/libifcore.a /opt/apps/intel/11.1/lib/intel64/libsvml.a")
else()
  set(MATH_LIBS "-L/opt/apps/intel/11.1/mkl/lib/em64t -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lguide -lpthread /opt/apps/intel/11.1/lib/intel64/libifcore.a /opt/apps/intel/11.1/lib/intel64/libsvml.a")
endif()
