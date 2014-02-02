### `src/lapack-like/HermitianTridiag/`

This folder contains Elemental's code for reducing a Hermitian matrix to 
real symmetric tridiagonal form. The various pieces are organized as follows:

-  `L.hpp`: Lower-triangular storage
-  `LSquare.hpp`: Lower-triangular storage specialized to
   square process grids
-  `LPan.hpp`: Panel portion of a blocked algorithm for lower-triangular 
   storage
-  `LPanSquare.hpp`: Panel portion of a blocked algorithm for lower-triangular
   storage specialized to square process grids
-  `U.hpp`: Upper-triangular storage
-  `USquare.hpp`: Upper-triangular storage specialized to square process grids
-  `UPan.hpp`: Panel portion of a blocked algorithm for upper-triangular 
   storage
-  `UPanSquare.hpp`: Panel portion of a blocked algorithm for upper-triangular
   storage specialized to square process grids
