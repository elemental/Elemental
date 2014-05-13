### `include/El/core/imports/`

This folder contains the header files for Elemental's wrappers
for external libraries, such as BLAS, LAPACK, and PMRRR. Please see 
`src/core/imports` for the corresponding source files.

In addition to this file, the folder contains:

-  `blas.hpp`: wrappers for Basic Linear Algebra Subprograms (BLAS)
-  `choice.hpp`: command-line processing for sequential programs
-  `flame.hpp`: wrappers for FLAME's QR-based bidiagonal SVD
-  `lapack.hpp`: wrappers for Linear Algebra PACKage (LAPACK)
-  `mpi.hpp`: wrappers for the Message Passing Interface (MPI)
-  `mpi_choice.hpp`: command-line processing for MPI-based programs
-  `pmrrr.hpp`: wrappers for Parallel Multiple Relatively Robust Representations
   (PMRRR)
