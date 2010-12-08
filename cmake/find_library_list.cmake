# A simple extension of find_library from a single lib to a list of libs
#
# Due to the find_library routine caching the return status for fixed variable 
# names within a loop, we must change it for each library.
macro(find_library_list FOUND_LIST SEARCH_LIST PATH PATH_SUFFIX)
  set(${FOUND_LIST})
  foreach(LIB ${${SEARCH_LIST}})
    find_library(FOUND_${LIB} ${LIB} PATHS ${PATH} PATH_SUFFIXES ${PATH_SUFFIX})
    if(FOUND_${LIB})
      list(APPEND ${FOUND_LIST} ${FOUND_${LIB}})
    endif(FOUND_${LIB})
    # On some machines, all of the above FOUND_${LIB} variables show up in 
    # the CMake wizard's cache. This is annoying, so we should explicitly 
    # remove them by forcing them to be internal.
    set(FOUND_${LIB} CACHE INTERNAL "")
  endforeach(LIB)
  # If FOUND_LIST is not the same length as SEARCH_LIST, then empty it so that
  # it is easy to check whether or not all of the libraries were found.
  list(LENGTH ${SEARCH_LIST} NUM_SEARCHED)
  list(LENGTH ${FOUND_LIST} NUM_FOUND)
  if(NUM_FOUND LESS NUM_SEARCHED)
    set(${FOUND_LIST})
  endif(NUM_FOUND LESS NUM_SEARCHED)
endmacro(find_library_list)

