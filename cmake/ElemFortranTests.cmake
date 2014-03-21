# Go ahead and check for Fortran, but keep in mind that CMake's 'OPTIONAL' 
# argument for enable_language is still completely broken as of 2.8.8
workaround_9220(Fortran FORTRAN_WORKS)
if(FORTRAN_WORKS)
  enable_language(Fortran OPTIONAL)
  set(HAVE_F90_INTERFACE FALSE)
  if(CMAKE_Fortran_COMPILER_SUPPORTS_F90)
    include(FortranCInterface)
    FortranCInterface_VERIFY(CXX)
    if(FortranCInterface_VERIFIED_CXX)
      set(HAVE_F90_INTERFACE TRUE)
      FortranCInterface_HEADER(
        ${CMAKE_CURRENT_BINARY_DIR}/include/elemental/FCMangle.h 
        MACRO_NAMESPACE "FC_")
      install(FILES ${PROJECT_BINARY_DIR}/include/elemental/FCMangle.h
              DESTINATION include/elemental/)
    endif()
  else()
    message(STATUS "${CMAKE_Fortran_COMPILER} does not appear to support F90")
  endif()
else()
  message(STATUS "Could not find working Fortran compiler")
endif()
