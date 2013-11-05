### `src/core/imports/`

This folder contains the directly-instantiated portions of Elemental's wrappers
for external libraries, such as BLAS, LAPACK, and PMRRR. Please see 
`include/elemental/core/imports` for the corresponding header-level prototypes
and implementations.

In addition to this file, the folder contains:

-  `blas.cpp`: wrappers for Basic Linear Algebra Subprograms (BLAS)
-  `flame.cpp`: wrappers for FLAME's QR-based bidiagonal SVD
-  `lapack.cpp`: wrappers for Linear Algebra PACKage (LAPACK)
-  `mpi.cpp`: wrappers for the Message Passing Interface (MPI)
-  `pmrrr.cpp`: wrappers for Parallel Multiple Relatively Robust Representations
   (PMRRR)
