### Current directory

This is the root of Elemental's `src/blas-like/Trr2k` folder, which contains 
the underlying implementations of rank-2k triangular updates.
In particular, in addition to this file, it holds:

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
