# The serial GNU compilers
set(COMPILER_DIR /opt/apps/gcc/4.7.1/bin)
set(CMAKE_C_COMPILER       ${COMPILER_DIR}/gcc)
set(CMAKE_CXX_COMPILER     ${COMPILER_DIR}/g++)
set(CMAKE_Fortran_COMPILER ${COMPILER_DIR}/gfortran)

# The MPI wrappers for the C and C++ compilers
set(MPI_COMPILER_DIR /opt/apps/gcc4_7/mvapich2/1.9/bin)
set(MPI_C_COMPILER       ${MPI_COMPILER_DIR}/mpicc)
set(MPI_CXX_COMPILER     ${MPI_COMPILER_DIR}/mpicxx)
set(MPI_Fortran_COMPILER ${MPI_COMPILER_DIR}/mpif90)

#set(MPI_C_COMPILE_FLAGS "")
#set(MPI_CXX_COMPILE_FLAGS "")
#set(MPI_Fortran_COMPILE_FLAGS "")
#set(MPI_C_INCLUDE_PATH       /opt/apps/gcc4_7/mvapich2/1.9/include)
#set(MPI_CXX_INCLUDE_PATH     ${MPI_C_INCLUDE_PATH})
#set(MPI_Fortran_INCLUDE_PATH ${MPI_C_INCLUDE_PATH})
#set(MPI_C_LINK_FLAGS "-Wl,-rpath,/opt/apps/limic2/0.5.5/lib -L/opt/apps/limic2/0.5.5/lib -L/opt/apps/gcc4_7/mvapich2/1.9/lib -L/opt/ofed/lib64/")
#set(MPI_CXX_LINK_FLAGS ${MPI_C_LINK_FLAGS})
#set(MPI_Fortran_LINK_FLAGS ${MPI_C_LINK_FLAGS})
#set(MPI_BASE_LIBS 
#    "-lmpich -lopa -llimic2 -lpthread -lrdmacm -libverbs -libumad -ldl -lrt")
#set(MPI_C_LIBRARIES "${MPI_BASE_LIBS}")
#set(MPI_CXX_LIBRARIES "-lmpichcxx ${MPI_BASE_LIBS}")
#set(MPI_Fortran_LIBRARIES "-lmpichf90 ${MPI_BASE_LIBS}")

if(CMAKE_BUILD_TYPE MATCHES PureDebug OR
   CMAKE_BUILD_TYPE MATCHES HybridDebug)
  set(CXX_FLAGS "-g -std=c++11")
else()
  set(CXX_FLAGS "-O3 -std=c++11")
endif()

set(OpenMP_CXX_FLAGS "-fopenmp")

set(MATH_LIBS "-L/opt/apps/intel/13/composer_xe_2013.2.146/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_lapack95_lp64 -lmkl_blas95_lp64 -lmkl_sequential -lpthread /opt/apps/gcc/4.7.1/lib64/libgfortran.a /opt/apps/intel/13/composer_xe_2013.2.146/compiler/lib/intel64/libimf.a /opt/apps/intel/13/composer_xe_2013.2.146/compiler/lib/intel64/libirc.a -lm")
