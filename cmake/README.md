### `cmake/`

This folder contains Elemental's auxiliary files for its CMake configuration
system:

-  `config.h.cmake`: The preprocessor definitions passed to Elemental
-  `ElemSub.cmake`: A CMake include file for simplifying the usage of
   Elemental as a CMake subproject
-  `ElemVars.cmake`: A Makefile include meant to simplify the usage of
   Elemental as a dependency in a project directly using 'make'
-  `tests/`: Various configuration tests
-  `toolchains/`: Toolchain files for various specialized architectures 
   (e.g., Blue Gene/Q)
