# Find Valgrind.
#
# This module defines:
#  VALGRIND_INCLUDE_DIR, where to find valgrind/memcheck.h, etc.
#  VALGRIND_PROGRAM, the valgrind executable.
#  VALGRIND_FOUND, If false, do not try to use valgrind.
#
# If you have valgrind installed in a non-standard place, you can define
# VALGRIND_PREFIX to tell cmake where it is.
#
# NOTE: Copied from the opencog project, where it is distributed under the
#       terms of the New BSD License.

message(STATUS "Valgrind Prefix: ${VALGRIND_PREFIX}")

find_path(VALGRIND_INCLUDE_DIR valgrind/memcheck.h
  /usr/include /usr/local/include ${VALGRIND_PREFIX}/include)
find_program(VALGRIND_PROGRAM NAMES valgrind PATH 
  /usr/bin /usr/local/bin ${VALGRIND_PREFIX}/bin)

find_package_handle_standard_args(VALGRIND DEFAULT_MSG
    VALGRIND_INCLUDE_DIR
    VALGRIND_PROGRAM)

mark_as_advanced(VALGRIND_INCLUDE_DIR VALGRIND_PROGRAM)
