### `src/blas_like/Trr2k/`

This folder contains the underlying implementations of rank-2k triangular 
updates. In particular, in addition to this file, it holds:

-  `Local.hpp`: all sequential implementations
-  `NNNN.hpp`: parallel normal/normal/normal/normal implementations
-  `NNNT.hpp`: parallel normal/normal/normal/transpose implementations
-  `NNTN.hpp`: parallel normal/normal/transpose/normal implementations
-  `NNTT.hpp`: etc.
-  `NTNN.hpp`
-  `NTNT.hpp`
-  `NTTN.hpp`
-  `NTTT.hpp`
-  `TNNN.hpp`
-  `TNNT.hpp`
-  `TNTN.hpp`
-  `TNTT.hpp`
-  `TTNN.hpp`
-  `TTNT.hpp`
-  `TTTN.hpp`
-  `TTTT.hpp`

#### Notes

There are currently two TODO items related to this folder:

1. Making the orientation options of LocalTrr2k more consistent with Trr2k
2. Implementing sequential versions of Trr2k
