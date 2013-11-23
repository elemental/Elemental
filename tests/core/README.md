### `tests/core`

This folder stores the correctness tests for Elemental's core functionality:

-  `AxpyInterface.cpp`: Tests the local-to-global and global-to-local Axpy 
   (y := alpha x plus y)  interface
-  `DifferentGrids.cpp`: Tests a redistribution between different process grids
-  `DistMatrix.cpp`: Tests various redistributions for the DistMatrix class
-  `Matrix.cpp`: Tests buffer attachment for the Matrix class
-  `Version.cpp`: Prints the version information of this Elemental build
