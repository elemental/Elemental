### `include/elemental/convex/`

- `decl.hpp`: include all declarations of convex opt. related routines
- `impl.hpp`: include all implementations of convex opt. related routines

This folder contains a few utilities for convex optimization; the only 
nontrivial one is for Singular Value soft-Thresholding (SVT).

-  `LogBarrier.hpp`: negative log of the determinant of an HPD matrix
-  `LogDetDiv.hpp`: divergence between two HPD matrices
-  `SoftThreshold.hpp`: soft-threshold each entry of a matrix
-  `SVT.hpp`: soft-threshold the singular values of a matrix
