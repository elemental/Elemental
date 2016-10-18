
function(MakeString ITEMS ITEMSTRING)
  set(${ITEMSTRING})
  foreach(ITEM ${${ITEMS}})
    set(${ITEMSTRING} "${${ITEMSTRING}} ${ITEM}")
  endforeach()
  set(${ITEMSTRING} ${${ITEMSTRING}} PARENT_SCOPE)
endfunction()

function(MakeLibString LIBS LIBSTRING)
  MakeString(${LIBS} ${LIBSTRING})
  set(${LIBSTRING} ${${LIBSTRING}} PARENT_SCOPE)
endfunction()

function(MakeIncString INCPATHS INCSTRING)
  set(${INCSTRING} PARENT_SCOPE)
  foreach(INCPATH ${${INCPATHS}})
    set(${INCSTRING} "${${INCSTRING}} -I${INCPATH}")
  endforeach()
  set(${INCSTRING} ${${INCSTRING}} PARENT_SCOPE)
endfunction()

MakeLibString(MATH_LIBS MATH_LIBSTRING)
MakeIncString(MPI_CXX_INCLUDE_PATH MPI_CXX_INCSTRING)
MakeLibString(MPI_CXX_LIBRARIES MPI_CXX_LIBSTRING)
MakeLibString(EXTERNAL_LIBS EXTERNAL_LIBSTRING)

MakeIncString(QD_INCLUDES QD_INCSTRING)
MakeIncString(MPC_INCLUDES MPC_INCSTRING)
MakeIncString(MPFR_INCLUDES MPFR_INCSTRING)
MakeIncString(GMP_INCLUDES GMP_INCSTRING)

MakeString(Qt5Widgets_DEFINITIONS QT5_DEFSTRING)
MakeIncString(Qt5Widgets_INCLUDE_DIRS QT5_INCSTRING)
MakeString(Qt5Widgets_EXECUTABLE_COMPILE_FLAGS QT5_COMPILESTRING)

# TODO: Generalize this for non-Unix architectures
set(QT5_LIBSTRING "-L${Qt5_LIBDIR} -lQt5Widgets -lQt5Gui -lQt5Core")
configure_file(${PROJECT_SOURCE_DIR}/cmake/configure_files/ElVars.in
               ${PROJECT_BINARY_DIR}/conf/ElVars @ONLY)
install(FILES ${PROJECT_BINARY_DIR}/conf/ElVars DESTINATION conf)
