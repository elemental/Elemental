### `src/lapack-like/`

This folder contains Elemental's source-level implementations of LAPACK-like
routines. Most such routines are implemented in `include/elemental/lapack-like`,
but the exceptions are `HermitianEig` and `HermitianTridiag`, which both involve
a substantial amount of code.

In addition to this file, this folder contains:

-  `HermitianEig.cpp`: Implementation of Hermitian eigensolvers
-  `HermitianTridiag/`: Underlying implementations of the reduction of a 
   Hermitian matrix to real symmetric tridiagonal form
-  `HermitianTridiag.cpp`: High-level interface for the reduction of a Hermitian
   matrix to real symmetric tridiagonal form
