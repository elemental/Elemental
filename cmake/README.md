### `cmake/`

This folder contains Elemental's auxiliary files for its CMake configuration
system:

-  `config.h.cmake`: The preprocessor definitions passed to Elemental
-  `elemental_sub.cmake`: A CMake include file for simplifying the usage of
   Elemental as a CMake subproject
-  `elemvariables.cmake`: A Makefile include meant to simplify the usage of
   Elemental as a dependency in a project directly using 'make'
-  `FindNumPy.cmake`: A CMake script for searching for the NumPy Python library
-  `GetGitRevisionDescription.cmake`: Main script for retrieving the current
   Git revision of Elemental
-  `GetGitRevisionDescription.cmake.in`: Input file as part of the script for
   determining the Git revision of Elemental
-  `language_support_v2.cmake`: A workaround for a CMake bug with respect to 
   checking for optional Fortran support
-  `toolchains/`: Toolchain files for various specialized architectures 
   (e.g., Blue Gene/Q)
