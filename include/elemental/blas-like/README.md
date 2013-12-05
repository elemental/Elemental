### `include/elemental/blas-like/`

- `decl.hpp`: include all declarations of BLAS-like routines
- `impl.hpp`: include all implementations of BLAS-like routines

Elemental's BLAS-like functionality is categorized into the typical 'levels':

-  level 1: operations with essentially no data reuse (e.g., DOT)
-  level 2: matrix/vector-like operations (e.g., GEMV)
-  level 3: matrix/matrix-like operations (e.g., GEMM)
