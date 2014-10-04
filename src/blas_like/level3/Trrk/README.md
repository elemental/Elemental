### `src/blas_like/Trrk/`

This folder contains the underlying implementations of rank-k triangular 
updates. In particular, in addition to this file, it holds:

-  `Local.hpp`: all sequential implementations
-  `NN.hpp`: parallel normal/normal implementations
-  `NT.hpp`: parallel normal/transposed implementations
-  `TN.hpp`: parallel transposed/normal implementations
-  `TT.hpp`: parallel transposed/transposed implementations

#### Notes

There is currently one TODO item related to this folder:

1. Making the orientation options of LocalTrrk more consistent with Trrk
