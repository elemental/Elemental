set(HAVE_QT5 FALSE)
if(Qt5_DIR)
  set(Qt5_LIBDIR "${Qt5_DIR}/lib")
endif()
if(Qt5_LIBDIR)
  set(Qt5Widgets_DIR "${Qt5_LIBDIR}/cmake/Qt5Widgets")
endif()
if(USE_QT5 OR Qt5Widgets_DIR)
  # Search for Qt5
  find_package(Qt5Widgets)
  if(Qt5Widgets_FOUND)
    set(ELEM_HEADERS_PREMOC 
        "include/elemental/io/display_window-premoc.hpp;include/elemental/io/complex_display_window-premoc.hpp")
    qt_wrap_cpp(elemental ELEM_MOC_SRC ${ELEM_HEADERS_PREMOC})
    include_directories(${Qt5Widgets_INCLUDE_DIRS})
    add_definitions(${Qt5Widgets_DEFINITIONS})
    set(EXTRA_FLAGS "${Qt5Widgets_EXECUTABLE_COMPILE_FLAGS} ${EXTRA_FLAGS}")

    # Qt5Widgets_DIR = Qt5_LIBDIR/cmake/Qt5Widgets
    get_filename_component(Qt5_CMAKEDIR ${Qt5Widgets_DIR} PATH)
    get_filename_component(Qt5_LIBDIR ${Qt5_CMAKEDIR} PATH)

    set(HAVE_QT5 TRUE)
    message(STATUS "Found Qt5")
  else()
    message(STATUS "Did NOT find Qt5")
  endif()
endif()
