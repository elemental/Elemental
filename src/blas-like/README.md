### `src/blas-like/`

This folder contains Elemental's source files for BLAS-like routines. 
The vast majority are in the `include/elemental/blas-like` folder, but 
the following are directly instantiated:

-  `Trr2k/`: underlying implementations of rank-2k triangular updates
-  `Trr2k.cpp`: the high-level interface to rank-2k triangular updates
-  `Trrk/`: underlying implementations of rank-k triangular updates
-  `Trrk.cpp`: the high-level interface to rank-k triangular updates
