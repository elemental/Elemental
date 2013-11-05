### `include/elemental/`

The structure of the header files in this folder follows a particular pattern.
For instance, for BLAS-like functionality, the relevant header files are as
follows:

-  `blas-like/`: supporting header-files for BLAS-like functionality
-  `blas-like_decl.hpp`: declarations for BLAS-like functionality implemented
   within `src/`
-  `blas-like.hpp`: includes all necessary header files for BLAS-like 
   functionality
-  `blas-like_impl.hpp`: includes only the implementations of BLAS-like 
   functionality

This pattern is followed for the majority of the following categories:

-  `blas-like`: BLAS-like functionality (e.g., GEMM)
-  `control`: solvers for control theory (e.g., Sylvester equations)
-  `convex/`: convex optimization routines (e.g., SVT)
-  `core/`: core data structures (e.g., `DistMatrix`)
-  `io/`: input/output functionality (e.g., printing and visualization)
-  `lapack-like/`: LAPACK-like functionality (e.g., SVD)
-  `matrices/`: special matrices (e.g., uniform random, Wilkinson, etc.)
