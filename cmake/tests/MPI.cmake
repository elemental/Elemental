include(CheckFunctionExists)
include(CheckCSourceCompiles)

find_package(MPI)
if(NOT MPI_C_FOUND)
  message(FATAL_ERROR "MPI C compiler was not found and is required")
endif()
include_directories(${MPI_C_INCLUDE_PATH})
set(EXTRA_FLAGS "${EXTRA_FLAGS} ${MPI_C_COMPILE_FLAGS}")

# NOTE: 
# check_function_exists only supports cdecl calling conventions, despite the
# fact that MS-MPI uses __stdcall. Thus, it is best to avoid using 
# check_function_exists in favor of check_c_source_compiles.
# Thanks to Ahn Vo for discovering this issue!
#
# TODO: Propagate the usage of check_c_source_compiles beyond MPI_Reduce_scatter
#       to the rest of this file.

# Ensure that we have MPI1 by looking for MPI_Reduce_scatter
# ==========================================================
set(CMAKE_REQUIRED_FLAGS "${MPI_C_COMPILE_FLAGS} ${MPI_LINK_FLAGS}")
set(CMAKE_REQUIRED_INCLUDES ${MPI_C_INCLUDE_PATH})
set(CMAKE_REQUIRED_LIBRARIES ${MPI_C_LIBRARIES})
set(MPI_REDUCE_SCATTER_CODE
    "#include \"mpi.h\"
     int main( int argc, char* argv[] )
     {
         MPI_Init( &argc, &argv );
         int *recvCounts;
         double *a, *b;  
         MPI_Reduce_scatter
         ( a, b, recvCounts, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
         MPI_Finalize();
         return 0;
     }")
check_c_source_compiles("${MPI_REDUCE_SCATTER_CODE}" EL_HAVE_MPI_REDUCE_SCATTER)
if(NOT EL_HAVE_MPI_REDUCE_SCATTER)
  message(FATAL_ERROR "Could not find MPI_Reduce_scatter")
endif()


# Test for whether or not we have Fortran MPI support
# ===================================================
if(MPI_Fortran_FOUND)
  set(MPIF_EXISTS FALSE)
  foreach(FORTRAN_INCLUDE_PATH ${MPI_Fortran_INCLUDE_PATH})
    if(EXISTS "${MPI_Fortran_INCLUDE_PATH}/mpif.h")
      set(MPIF_EXISTS TRUE)
    endif()
  endforeach()
  if(MPIF_EXISTS)
    set(EL_HAVE_MPI_FORTRAN TRUE)
  else()
    message(WARNING 
      "Fortran MPI support detected, but mpif.h was not found in ${MPI_Fortran_INCLUDE_PATH}")
  endif()
endif()

# Ensure that we have MPI2 by looking for MPI_Type_create_struct
# ==============================================================
check_function_exists(MPI_Type_create_struct EL_HAVE_MPI_TYPE_CREATE_STRUCT)
if(NOT EL_HAVE_MPI_TYPE_CREATE_STRUCT)
  message(FATAL_ERROR "Could not find MPI_Type_create_struct")
endif()

# Ensure that we have MPI_[UNSIGNED_]LONG_LONG if compiling with 64-bit integers
# ==============================================================================
set(MPI_LONG_LONG_CODE
    "#include \"mpi.h\"
     int main( int argc, char* argv[] )
     {
         MPI_Init( &argc, &argv );
         MPI_Datatype lli = MPI_LONG_LONG_INT;
         MPI_Datatype llu = MPI_UNSIGNED_LONG_LONG;
         MPI_Finalize();
         return 0;
     }")
check_c_source_compiles("${MPI_LONG_LONG_CODE}" EL_HAVE_MPI_LONG_LONG)
if(EL_USE_64BIT_INTS AND NOT EL_HAVE_MPI_LONG_LONG)
  message(FATAL_ERROR 
    "Did not detect MPI_LONG_LONG_INT and MPI_UNSIGNED_LONG_LONG")
endif()

# Check if MPI_LONG_DOUBLE[_COMPLEX] is supported
# ===============================================
set(MPI_LONG_DOUBLE_CODE
    "#include \"mpi.h\"
     int main( int argc, char* argv[] )
     {
         MPI_Init( &argc, &argv );
         MPI_Datatype longDbl = MPI_LONG_DOUBLE;
         MPI_Finalize();
         return 0;
     }")
set(MPI_LONG_DOUBLE_COMPLEX_CODE
    "#include \"mpi.h\"
     int main( int argc, char* argv[] )
     {
         MPI_Init( &argc, &argv );
         MPI_Datatype longDblCpx = MPI_LONG_DOUBLE_COMPLEX;
         MPI_Finalize();
         return 0;
     }")
check_c_source_compiles("${MPI_LONG_DOUBLE_CODE}" EL_HAVE_MPI_LONG_DOUBLE)
check_c_source_compiles("${MPI_LONG_DOUBLE_COMPLEX_CODE}" 
  EL_HAVE_MPI_LONG_DOUBLE_COMPLEX)

# Check if MPI_C_FLOAT_COMPLEX and MPI_C_DOUBLE_COMPLEX exist
# ===========================================================
set(MPI_C_COMPLEX_CODE
    "#include \"mpi.h\"
     int main( int argc, char* argv[] )
     {
         MPI_Init( &argc, &argv );
         MPI_Datatype floatCpx = MPI_C_FLOAT_COMPLEX;
         MPI_Datatype doubleCpx = MPI_C_DOUBLE_COMPLEX;
         MPI_Finalize();
         return 0;
     }")
check_c_source_compiles("${MPI_C_COMPLEX_CODE}" EL_HAVE_MPI_C_COMPLEX)

# Detect support for various optional MPI routines
# ================================================
check_function_exists(MPI_Reduce_scatter_block EL_HAVE_MPI_REDUCE_SCATTER_BLOCK)
check_function_exists(MPI_Iallgather EL_HAVE_MPI3_NONBLOCKING_COLLECTIVES)
check_function_exists(MPIX_Iallgather EL_HAVE_MPIX_NONBLOCKING_COLLECTIVES)
check_function_exists(MPI_Init_thread EL_HAVE_MPI_INIT_THREAD)
check_function_exists(MPI_Query_thread EL_HAVE_MPI_QUERY_THREAD)
check_function_exists(MPI_Comm_set_errhandler EL_HAVE_MPI_COMM_SET_ERRHANDLER)
# Detecting MPI_IN_PLACE and MPI_Comm_f2c requires test compilation
# -----------------------------------------------------------------
set(MPI_IN_PLACE_CODE
    "#include \"mpi.h\"
     int main( int argc, char* argv[] )
     {
         MPI_Init( &argc, &argv );
         float a;
         MPI_Allreduce
         ( MPI_IN_PLACE, &a, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD );
         MPI_Finalize();
         return 0;
     }")
check_c_source_compiles("${MPI_IN_PLACE_CODE}" EL_HAVE_MPI_IN_PLACE)
set(MPI_COMM_F2C_CODE
"#include \"mpi.h\"
 int main( int argc, char* argv[] )
 {
     MPI_Init( &argc, &argv );
     MPI_Fint fint;
     MPI_Comm comm = MPI_Comm_f2c( fint );
     MPI_Finalize();
     return 0;
 }")
check_c_source_compiles("${MPI_COMM_F2C_CODE}" EL_HAVE_MPI_COMM_F2C)

# Detect whether or not MPI_Comm and MPI_Group are implemented as an int
# ======================================================================
# NOTE: These are no longer used internally by Elemental 
set(MPI_COMM_NOT_INT_CODE
    "#include \"mpi.h\"
     void Foo( MPI_Comm comm ) { }
     void Foo( int comm ) { }
     int main( int argc, char* argv[] )
     {
         MPI_Init( &argc, &argv );
         MPI_Finalize();
         return 0;
     }")
set(MPI_GROUP_NOT_INT_CODE
    "#include \"mpi.h\"
     void Foo( MPI_Group group ) { }
     void Foo( int group ) { }
     int main( int argc, char* argv[] )
     {
         MPI_Init( &argc, &argv );
         MPI_Finalize();
         return 0;
     }")
check_c_source_compiles("${MPI_COMM_NOT_INT_CODE}" EL_MPI_COMM_NOT_INT)
check_c_source_compiles("${MPI_GROUP_NOT_INT_CODE}" EL_MPI_GROUP_NOT_INT)
