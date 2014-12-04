find_package(MPI)
if(NOT MPI_CXX_FOUND)
  message(FATAL_ERROR "MPI C++ compiler was not found and is required")
endif()
include_directories(${MPI_CXX_INCLUDE_PATH})
set(EXTRA_FLAGS "${EXTRA_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
set(CMAKE_REQUIRED_FLAGS "${MPI_CXX_COMPILE_FLAGS} ${MPI_LINK_FLAGS}")
set(CMAKE_REQUIRED_INCLUDES ${MPI_CXX_INCLUDE_PATH})
set(CMAKE_REQUIRED_LIBRARIES ${MPI_CXX_LIBRARIES})

# Ensure that we have MPI1 by looking for MPI_Reduce_scatter
# ==========================================================
check_function_exists(MPI_Reduce_scatter EL_HAVE_MPI_REDUCE_SCATTER)
if(NOT EL_HAVE_MPI_REDUCE_SCATTER)
  message(FATAL_ERROR "Could not find MPI_Reduce_scatter")
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
check_cxx_source_compiles("${MPI_LONG_LONG_CODE}" EL_HAVE_MPI_LONG_LONG)
if(EL_USE_64BIT_INTS AND NOT EL_HAVE_MPI_LONG_LONG)
  message(FATAL_ERROR 
    "Did not detect MPI_LONG_LONG_INT and MPI_UNSIGNED_LONG_LONG")
endif()

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
check_cxx_source_compiles("${MPI_IN_PLACE_CODE}" EL_HAVE_MPI_IN_PLACE)
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
check_cxx_source_compiles("${MPI_COMM_F2C_CODE}" EL_HAVE_MPI_COMM_F2C)

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
check_cxx_source_compiles("${MPI_COMM_NOT_INT_CODE}" EL_MPI_COMM_NOT_INT)
check_cxx_source_compiles("${MPI_GROUP_NOT_INT_CODE}" EL_MPI_GROUP_NOT_INT)
