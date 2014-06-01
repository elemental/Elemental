# MPI wrappers assumed to be in path
set(COMPILER_DIR )
set(CMAKE_C_COMPILER       ${COMPILER_DIR}mpicc)
set(CMAKE_CXX_COMPILER     ${COMPILER_DIR}mpicxx)
set(CMAKE_Fortran_COMPILER ${COMPILER_DIR}mpif90)

set(MPI_C_COMPILER ${CMAKE_C_COMPILER})
set(MPI_CXX_COMPILER ${CMAKE_CXX_COMPILER})
set(MPI_Fortran_COMPILER ${CMAKE_Fortran_COMPILER})

if(CMAKE_BUILD_TYPE MATCHES PureDebug OR
   CMAKE_BUILD_TYPE MATCHES HybridDebug)
  set(C_FLAGS   "-g")
  set(CXX_FLAGS "-g")
else()
  set(C_FLAGS   "-O3")
  set(CXX_FLAGS "-O3")
endif()

set(OpenMP_CXX_FLAGS "-fopenmp")

set(MATH_LIBS "-framework Accelerate")
set(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
set(BUILD_SHARED_LIBS ON)
set(CMAKE_EXE_LINKER_FLAGS "")
