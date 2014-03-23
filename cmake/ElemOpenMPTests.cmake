if(HYBRID)
  if(OpenMP_CXX_FLAGS)
    set(ELEM_HAVE_OPENMP TRUE)
    set(OpenMP_CXX_FLAGS ${OpenMP_CXX_FLAGS})
    message(STATUS "Using OpenMP_CXX_FLAGS=${OpenMP_CXX_FLAGS}")
  else()
    find_package(OpenMP)
    if(OPENMP_FOUND)
      set(ELEM_HAVE_OPENMP TRUE)
    else()
      set(OpenMP_CXX_FLAGS "" CACHE STRING "OpenMP CXX FLAGS")
      message(FATAL_ERROR
        "Hybrid build failed because OpenMP support not detected. Please specify OpenMP_CXX_FLAGS.")
    endif()
  endif()
  set(EXTRA_FLAGS "${OpenMP_CXX_FLAGS} ${EXTRA_FLAGS}")
endif()
# See if we have 'collapse' support, since Clang is currently setting _OPENMP
# greater than 200805 despite not supporting OpenMP 3.0 (i.e., collapse)
set(CMAKE_REQUIRED_FLAGS ${OpenMP_CXX_FLAGS})
set(OMP_COLLAPSE_CODE
    "#include <omp.h>
     int main( int argc, char* argv[] ) 
     {
         int k[100];
     #pragma omp collapse(2)
         for( int i=0; i<10; ++i ) 
             for( int j=0; j<10; ++j )
                 k[i+j*10] = i+j; 
         return 0; 
     }")
check_cxx_source_compiles("${OMP_COLLAPSE_CODE}" ELEM_HAVE_OMP_COLLAPSE)
