### `src/`

This folder contains Elemental's routines which are to be directly instantiated 
rather than keeping them datatype-agnostic. Several commonly-used large routines
are kept here in order to keep the build times reasonable. This folder contains 
the following subfolders:

-  `blas-like/`: BLAS-like routines
-  `core/`: Elemental's core data structures and environment
-  `io/`: input/output, such as Qt5 graphics
-  `lapack-like/`: LAPACK-like routines
