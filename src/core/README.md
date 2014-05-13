### `src/core/`

This folder contains the directly-instantiated portions of Elemental's core
functionality, such as the `Matrix` and `DistMatrix` classes, and wrappers of
external libraries, such as BLAS and LAPACK. Please see 
`include/El/core` for the corresponding header-level implementations and
prototypes.

In addition to this file, the folder currently contains:

-  `dist_matrix/`: distributed matrix class (`DistMatrix`)
-  `global.cpp`: initialization/finalization, call-stack manipulation, etc.
-  `imports/`: wrappers for external software
-  `matrix.cpp`: sequential matrix class (`Matrix`)
-  `mpi_register.cpp`: custom MPI datatypes and reduction operations
