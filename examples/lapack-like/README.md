### `examples/lapack-like`

This folder contains several examples of Elemental's LAPACK-like functionality:

-  `BunchKaufman.cpp`: Accurate symmetric/Hermitian-indefinite factorization
-  `BusingerGolub.cpp`: Column-pivoted QR decomposition
-  `ComplexHermitianFunction.cpp`: Applies a complex function to the eigenvalues
   of a Hermitian matrix
-  `GaussianElimination.cpp`: Solves systems of equations via Gaussian elim.
-  `HermitianEig.cpp`: Computes the eigen{values/pairs} of a Hermitian matrix
-  `HermitianEigFromSequential.cpp`: Distributes a sequential Hermitian matrix, computes its EVD, and then gathers the result back to the original process
-  `HermitianPseudoinverse.cpp`: Forms the pseudoinverse of a Hermitian matrix
-  `HermitianQDWH.cpp`: A variant of the QDWH algorithm for the polar 
   decomposition which is specialized for Hermitian matrices
-  `HermitianSDC.cpp`: Spectral Divide and Conquer eigensolver for Hermitian 
   matrices
-  `HermitianSVD.cpp`: Singular Value Decomposition of a Hermitian matrix
-  `HPDInverse.cpp`: Inverts a Hermitian Positive-Definite matrix
-  `HPSDCholesky.cpp`: Computes the (non-unique) Cholesky decomposition of a 
   Hermitian Positive-SemiDefinite matrix via its eigenvalue decomposition
-  `HPSDSquareRoot.cpp`: Computes the square-root of a Hermitian 
   Positive-SemiDefinite matrix via its eigenvalue decomposition
-  `ID.cpp`: Computes an Interpolate Decomposition 
   (closely related to pivoted QR)
-  `KyFanAndSchatten.cpp`: Compute Ky Fan and Schatten norms
-  `LDL.cpp`: Unpivoted LDL^T/LDL^H factorization
-  `LDLInverse.cpp`: Invert a symmetric/Hermitian matrix via a pivoted 
   symmetric factorization (e.g., Bunch-Kaufman)
-  `LeastSquares.cpp`: Solve a least-squares problem via a QR decomposition
-  `Polar.cpp`: Compute a polar decomposition (unitary times HPD)
-  `Pseudoinverse.cpp`: Compute the pseudoinverse of an arbitrary matrix
-  `QDWH.cpp`: Compute the polar factor of an arbitrary matrix via the QDWH 
   algorithm
-  `QR.cpp`: Compute a QR decomposition
-  `RealHermitianFunction.cpp`: Apply a real function to the eigenvalues of a
   Hermitian matrix
-  `RealSchur.cpp`: Compute the Schur decomposition of a real matrix
-  `RealSymmetricFunction.cpp`: Apply a real function to the eigenvalues of a 
   (real) symmetric matrix
-  `Schur.cpp`: Compute the Schur decomposition of a matrix
-  `SequentialBunchKaufman.cpp`: Test the sequential algorithm for Bunch-Kaufman
-  `SequentialQR.cpp`: Test the sequential algorithm for QR decomposition
-  `SequentialSVD.cpp`: Test the sequential algorithm for SVD
-  `Sign.cpp`: Test the matrix sign function (maps eigenvalues to {-1,+1})
-  `SimpleSVD.cpp`: An extremely simple SVD driver
-  `Skeleton.cpp`: Compute a matrix skeleton
-  `SkewHermitianEig.cpp`: Compute the EVD of a skew-Hermitian matrix
-  `SVD.cpp`: Compute the SVD of an arbitrary matrix
