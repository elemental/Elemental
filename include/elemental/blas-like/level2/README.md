### `include/elemental/blas-like/level2/`

The level-2 BLAS-like routines implemented in this folder are:

-  `ApplyColumnPivots.hpp`: apply a permutation matrix from the right
-  `ApplyRowPivots.hpp`: apply a permutation matrix from the left
-  `ApplySymmetricPivots.hpp`: apply a symmetric permutation to a symmetric
   matrix
-  `ComposePivots.hpp`: form the explicit image and preimage of a permutation
-  `Gemv/`: supporting routines for GEneral Matrix/Vector multiplication (GEMV)
-  `Gemv.hpp`: interface for GEMV
-  `Ger.hpp`: GEneral Rank-one update (GER)
-  `Geru.hpp`: unconjugated version of GER
-  `Hemv.hpp`: HErmitian Matrix/Vector multiplication (HEMV)
-  `Her2.hpp`: HErmitian Rank-2 update (HER2)
-  `Her.hpp`: HErmitian Rank-1 update (HER)
-  `Symv/`: supporting routines for SYmmetric Matrix/Vector multiplication 
   (SYMV)
-  `Symv.hpp`: interface for SYMV
-  `Syr2.hpp`: SYmmetric Rank-2 update (SYR2)
-  `Syr.hpp`: SYmmetric Rank-1 update (SYR)
-  `Trmv.hpp`: TRiangular Matrix/Vector multiplication (TRMV)
-  `Trr2.hpp`: TRiangular Rank-2 update (TRR2)
-  `Trr.hpp`: TRiangular Rank-1 update (TRR)
-  `Trsv/`: supporting routines for TRiangular Solve against Vector (TRSV)
-  `Trsv.hpp`: interface for TRSV

NOTE: The symmetric pivot applications are not yet optimized.
